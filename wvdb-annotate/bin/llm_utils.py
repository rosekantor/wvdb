#!/usr/bin/env python3
"""
llm_utils.py — shared multi-vendor LLM call dispatcher.

Supports three provider modes via raw HTTP requests (requests.post), not
the vendor SDKs — this keeps the conda environment lean (no anthropic/
openai packages to version-pin) and makes the exact request/response
shape fully visible and easy to debug or extend to a new provider:

  anthropic          — public Anthropic API (api.anthropic.com)
  openai             — public OpenAI API (api.openai.com)
  openai_compatible  — any endpoint that implements the same
                       chat-completions request/response shape as OpenAI
                       (Authorization: Bearer, "messages" array in,
                       choices[0].message.content out). This covers
                       internal enterprise LLM gateways (a common pattern
                       — many organizations front multiple vendor models
                       behind one OpenAI-shaped endpoint), Azure OpenAI,
                       and locally-hosted servers (vLLM, Ollama, etc).
                       Requires --llm-base-url; there is no built-in
                       default URL for this mode.

Both provider and model must be specified explicitly by the caller —
there is no default model, since model availability and quality change
over time and a silently-used default could go stale without anyone
noticing.

API keys are loaded from a .env file (via python-dotenv), never from a
hardcoded path or committed default. The .env file's path is provided
explicitly by the calling script/pipeline, read directly from the
filesystem — it is never staged as a Nextflow `path` input, so it never
enters `work/` or any trace/log directory Nextflow manages.

Expected .env keys (only the one matching --provider is required):
  ANTHROPIC_API_KEY=sk-ant-...
  OPENAI_API_KEY=sk-...
  OPENAI_COMPATIBLE_API_KEY=...       (for an internal/custom gateway)

If your organization's convention uses a different .env variable name,
override it with --llm-api-key-env-var instead of one of the defaults
above.
"""

import json
import sys
import time

import requests
from dotenv import dotenv_values


SUPPORTED_PROVIDERS = ('anthropic', 'openai', 'openai_compatible')

ANTHROPIC_URL = "https://api.anthropic.com/v1/messages"
ANTHROPIC_VERSION = "2023-06-01"
OPENAI_URL = "https://api.openai.com/v1/chat/completions"

DEFAULT_API_KEY_ENV_VAR = {
    'anthropic':         'ANTHROPIC_API_KEY',
    'openai':             'OPENAI_API_KEY',
    'openai_compatible':  'OPENAI_COMPATIBLE_API_KEY',
}


def load_api_key(env_file, provider, api_key_env_var=None):
    """
    Load the API key for the given provider from a .env file.
    Raises a clear error if the file or key is missing, rather than
    silently proceeding with no authentication.

    api_key_env_var overrides the default .env variable name looked up
    for this provider — use this if your organization's gateway stores
    the key under a different name than the defaults above.
    """
    if provider not in SUPPORTED_PROVIDERS:
        raise ValueError(
            f"Unsupported provider: {provider!r}. "
            f"Supported: {SUPPORTED_PROVIDERS}"
        )

    env_key_name = api_key_env_var or DEFAULT_API_KEY_ENV_VAR[provider]

    values = dotenv_values(env_file)
    api_key = values.get(env_key_name)
    if not api_key:
        raise ValueError(
            f"{env_key_name} not found in {env_file}. "
            f"Add a line: {env_key_name}=your-key-here"
        )
    return api_key


def call_llm(provider, model, api_key, system_prompt, user_prompt,
             base_url=None, max_retries=3, retry_delay=5, max_tokens=4096):
    """
    Call the specified provider/model with a system + user prompt.
    Returns the raw text response. Retries on transient failures rather
    than crashing immediately — LLM API calls are exactly the kind of
    network-dependent operation that benefits from the same
    retry-then-skip-gracefully pattern used for Entrez lookups elsewhere
    in this pipeline.

    base_url is required when provider == 'openai_compatible'; ignored
    otherwise (anthropic/openai always use their fixed public endpoints).
    """
    if provider not in SUPPORTED_PROVIDERS:
        raise ValueError(
            f"Unsupported provider: {provider!r}. "
            f"Supported: {SUPPORTED_PROVIDERS}"
        )
    if provider == 'openai_compatible' and not base_url:
        raise ValueError(
            "--llm-base-url is required when --provider openai_compatible"
        )

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            if provider == 'anthropic':
                return _call_anthropic(model, api_key, system_prompt,
                                       user_prompt, max_tokens)
            elif provider == 'openai':
                return _call_openai_shaped(OPENAI_URL, model, api_key,
                                           system_prompt, user_prompt, max_tokens)
            elif provider == 'openai_compatible':
                return _call_openai_shaped(base_url, model, api_key,
                                           system_prompt, user_prompt, max_tokens)
        except Exception as e:
            last_error = e
            print(f"  LLM call attempt {attempt}/{max_retries} failed: {e}",
                  file=sys.stderr)
            if attempt < max_retries:
                time.sleep(retry_delay)

    raise RuntimeError(
        f"LLM call failed after {max_retries} attempts: {last_error}"
    )


def _call_anthropic(model, api_key, system_prompt, user_prompt, max_tokens):
    headers = {
        "x-api-key": api_key,
        "anthropic-version": ANTHROPIC_VERSION,
        "content-type": "application/json",
    }
    payload = {
        "model": model,
        "max_tokens": max_tokens,
        "system": system_prompt,
        "messages": [{"role": "user", "content": user_prompt}],
    }
    response = requests.post(ANTHROPIC_URL, headers=headers, json=payload, timeout=120)
    response.raise_for_status()
    data = response.json()
    return data["content"][0]["text"]


def _call_openai_shaped(url, model, api_key, system_prompt, user_prompt, max_tokens):
    """
    Shared implementation for both 'openai' and 'openai_compatible' — the
    public OpenAI API is itself just one instance of the OpenAI-shaped
    chat-completions interface, so the same request/response handling
    covers the public API, internal gateways, Azure OpenAI, and locally
    hosted OpenAI-compatible servers alike. Only the URL differs.

    Response parsing uses .get() with defaults rather than direct
    indexing, so a malformed/unexpected response shape (e.g. a gateway
    that wraps or filters the response slightly differently) raises a
    clear empty-reply condition rather than a raw KeyError.
    """
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }
    payload = {
        "model": model,
        "max_tokens": max_tokens,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    }
    response = requests.post(url, headers=headers, json=payload, timeout=120)
    response.raise_for_status()
    data = response.json()
    content = data.get("choices", [{}])[0].get("message", {}).get("content", "")
    if not content:
        raise RuntimeError(
            f"Empty or unexpected response shape from {url}: {data!r}"
        )
    return content


def parse_json_response(text):
    """
    Extract a JSON object from an LLM response, stripping markdown code
    fences if present (models often wrap JSON in ```json ... ``` even
    when explicitly asked not to).
    """
    text = text.strip()
    if text.startswith('```'):
        lines = text.split('\n')
        # Drop the opening fence (```json or ```) and closing fence
        lines = lines[1:]
        if lines and lines[-1].strip() == '```':
            lines = lines[:-1]
        text = '\n'.join(lines)
    return json.loads(text)
