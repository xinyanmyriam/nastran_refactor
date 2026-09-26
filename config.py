"""
NASTRAN Refactoring Experiment Configuration.

API 配置：
- 国产模型（DeepSeek, Qwen, Kimi）通过 OpenRouter
- Claude/GPT 通过 API2D 中转站
"""

# =====================================================
# API Keys
# =====================================================

OPENROUTER_API_KEY = "sk-or-v1-REPLACE_WITH_YOUR_OPENROUTER_KEY"

API2D_KEY = "sk-REPLACE_WITH_YOUR_API2D_KEY"
API2D_URL = "https://api.gptsapi.net"

DEEPSEEK_KEY = "sk-REPLACE_WITH_YOUR_DEEPSEEK_KEY"
DEEPSEEK_URL = "https://api.deepseek.com"

QWEN_KEY = "sk-REPLACE_WITH_YOUR_QWEN_KEY"
QWEN_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"

# =====================================================
# 模型配置
# =====================================================

MODELS = {
    # --- 闭源 Frontier（通过 API2D）---
    "claude-sonnet-4": {
        "provider": "openai",
        "model": "claude-sonnet-4-6",
        "api_key": API2D_KEY,
        "base_url": API2D_URL + "/v1",
        "max_tokens": 16384,
    },
    "gpt-5": {
        "provider": "openai",
        "model": "gpt-5",
        "api_key": API2D_KEY,
        "base_url": API2D_URL + "/v1",
        "max_tokens": 8192,
    },

    # --- 国产开源 Frontier ---
    "deepseek-v4-pro": {
        "provider": "openai",
        "model": "deepseek-v4-pro",
        "api_key": DEEPSEEK_KEY,
        "base_url": DEEPSEEK_URL + "/v1",
        "max_tokens": 16384,
    },
    "deepseek-v4-flash": {
        "provider": "openai",
        "model": "deepseek-v4-flash",
        "api_key": DEEPSEEK_KEY,
        "base_url": DEEPSEEK_URL + "/v1",
        "max_tokens": 16384,
    },
    "qwen3": {
        "provider": "openai",
        "model": "qwen/qwen3-235b-a22b",
        "api_key": OPENROUTER_API_KEY,
        "base_url": "https://openrouter.ai/api/v1",
        "max_tokens": 8192,
    },
    # Alibaba Tongyi Qwen (official DashScope OpenAI-compatible endpoint).
    # `model` is the DashScope API identifier, not the console display name;
    # confirm the exact snapshot string in the DashScope console if calls 404.
    "qwen3.7-plus": {
        "provider": "openai",
        "model": "qwen-plus",
        "api_key": QWEN_KEY,
        "base_url": QWEN_URL,
        "max_tokens": 16384,
    },
    "kimi-k2.5": {
        "provider": "openai",
        "model": "moonshotai/kimi-k2.5",
        "api_key": OPENROUTER_API_KEY,
        "base_url": "https://openrouter.ai/api/v1",
        "max_tokens": 8192,
    },
    "deepseek-r1": {
        "provider": "openai",
        "model": "deepseek-reasoner",
        "api_key": DEEPSEEK_KEY,
        "base_url": DEEPSEEK_URL + "/v1",
        "max_tokens": 16384,
    },
}

# 默认使用的模型（国产模型走 OpenRouter，更稳定）
DEFAULT_MODEL = "deepseek-v4-flash"

# =====================================================
# NASTRAN 源码路径
# =====================================================

NASTRAN_SOURCE_DIR = "nastran/NASTRAN-95/mis"
NASTRAN_INPUT_DIR = "nastran/NASTRAN-95/inp"
NASTRAN_OUTPUT_DIR = "nastran/NASTRAN-95/demoout"

