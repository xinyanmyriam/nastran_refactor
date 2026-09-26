"""
LLM Client for NASTRAN Refactoring Experiment.

Provides a unified interface for calling LLMs through OpenRouter.
Compatible with the multi-agent framework's expected interface:
    client.generate(system_prompt, user_prompt) -> str
"""

import time
import json
from openai import OpenAI

from config import MODELS, DEFAULT_MODEL


class LLMClient:
    """
    LLM client that wraps OpenRouter API calls.
    
    Usage:
        from llm_client import LLMClient
        client = LLMClient()  # Uses DEFAULT_MODEL
        response = client.generate(system_prompt, user_prompt)
    """
    
    def __init__(self, model: str = None, temperature: float = 0.0, max_retries: int = 3, logger=None):
        """
        Args:
            model: Model name key from config.MODELS. Defaults to config.DEFAULT_MODEL.
            temperature: Sampling temperature (0 = deterministic).
            max_retries: Number of retry attempts on failure.
            logger: Optional ExperimentLogger for recording LLM calls.
        """
        self.model_name = model or DEFAULT_MODEL
        self.logger = logger
        
        if self.model_name not in MODELS:
            available = ", ".join(MODELS.keys())
            raise ValueError(
                f"Unknown model '{self.model_name}'. Available: {available}"
            )
        
        self.config = MODELS[self.model_name]
        self.temperature = temperature
        self.max_retries = max_retries
        self.max_tokens = self.config.get("max_tokens", 4096)
        
        # Initialize OpenAI client (OpenRouter is OpenAI-compatible)
        self._client = OpenAI(
            api_key=self.config["api_key"],
            base_url=self.config.get("base_url", "https://openrouter.ai/api/v1"),
            timeout=900.0,  # 15 minutes timeout for R1 reasoning model
        )
        
        # Statistics
        self.total_calls = 0
        self.total_tokens = 0
    
    def generate(self, system_prompt: str, user_prompt: str, agent_name: str = "unknown") -> str:
        """
        Generate a response from the LLM.
        
        This is the main interface used by all agents.
        
        Args:
            system_prompt: System message defining the agent's role.
            user_prompt: User message with the specific task.
            agent_name: Name of the calling agent (for logging).
            
        Returns:
            str: The model's response text.
        """
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ]
        
        call_start = time.time()
        
        # Build extra parameters for DeepSeek thinking mode control
        extra_kwargs = {}
        if "deepseek" in self.config.get("base_url", ""):
            # Disable thinking mode for code generation (avoids reasoning token exhaustion)
            extra_kwargs["extra_body"] = {"thinking": {"type": "disabled"}}
        
        for attempt in range(self.max_retries):
            try:
                response = self._client.chat.completions.create(
                    model=self.config["model"],
                    messages=messages,
                    temperature=self.temperature,
                    max_tokens=self.max_tokens,
                    **extra_kwargs,
                )
                
                self.total_calls += 1
                tokens_used = 0
                if response.usage:
                    tokens_used = response.usage.total_tokens
                    self.total_tokens += tokens_used
                
                content = response.choices[0].message.content
                result = content.strip() if content else ""
                
                elapsed = time.time() - call_start
                
                # Log the call
                if self.logger:
                    self.logger.log_llm_call(
                        agent_name=agent_name,
                        system_prompt=system_prompt,
                        user_prompt=user_prompt,
                        response=result,
                        elapsed_seconds=elapsed,
                        tokens_used=tokens_used,
                    )
                
                return result
                
            except Exception as e:
                elapsed = time.time() - call_start
                print(f"  [LLM] Attempt {attempt + 1}/{self.max_retries} failed: {e}")
                
                if self.logger:
                    self.logger.log_llm_call(
                        agent_name=agent_name,
                        system_prompt=system_prompt,
                        user_prompt=user_prompt,
                        response="",
                        elapsed_seconds=elapsed,
                        tokens_used=0,
                        error=str(e),
                    )
                
                if attempt < self.max_retries - 1:
                    wait = 2 ** attempt
                    print(f"  [LLM] Retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    # Final retry with longer wait
                    print(f"  [LLM] All {self.max_retries} retries exhausted. Waiting 30s for final attempt...")
                    time.sleep(30)
                    try:
                        response = self._client.chat.completions.create(
                            model=self.config["model"],
                            messages=messages,
                            temperature=self.temperature,
                            max_tokens=self.max_tokens,
                            **extra_kwargs,
                        )
                        content = response.choices[0].message.content
                        result = content.strip() if content else ""
                        if result and not result.startswith("ERROR"):
                            self.total_calls += 1
                            elapsed = time.time() - call_start
                            if self.logger:
                                self.logger.log_llm_call(
                                    agent_name=agent_name,
                                    system_prompt=system_prompt,
                                    user_prompt=user_prompt,
                                    response=result,
                                    elapsed_seconds=elapsed,
                                    tokens_used=response.usage.total_tokens if response.usage else 0,
                                )
                            return result
                    except Exception as e2:
                        pass
                    print(f"  [LLM] All retries exhausted.")
                    # Log with clear API_TIMEOUT marker for later identification
                    if self.logger:
                        self.logger.log_llm_call(
                            agent_name=agent_name,
                            system_prompt="",
                            user_prompt="",
                            response="",
                            elapsed_seconds=time.time() - call_start,
                            tokens_used=0,
                            error=f"API_TIMEOUT: All {self.max_retries + 1} attempts failed. Last error: {e}",
                        )
                    return ""
    
    def generate_simple(self, prompt: str) -> str:
        """
        Simple generation without system prompt.
        
        Args:
            prompt: The prompt text.
            
        Returns:
            str: The model's response.
        """
        return self.generate(
            system_prompt="You are a helpful assistant.",
            user_prompt=prompt
        )
    
    def get_stats(self) -> dict:
        """Get usage statistics."""
        return {
            "model": self.model_name,
            "total_calls": self.total_calls,
            "total_tokens": self.total_tokens,
        }
    
    def __repr__(self):
        return f"LLMClient(model='{self.model_name}', calls={self.total_calls})"


# Convenience function for quick testing
def quick_test():
    """Quick test to verify API connectivity."""
    client = LLMClient()
    print(f"Testing model: {client.model_name}")
    print(f"Config: {client.config['model']} via {client.config.get('base_url', 'default')}")
    
    response = client.generate(
        system_prompt="You are a helpful assistant.",
        user_prompt="Say 'API connection successful' and nothing else."
    )
    
    print(f"Response: {response}")
    print(f"Stats: {client.get_stats()}")
    return "successful" in response.lower()


if __name__ == "__main__":
    success = quick_test()
    if success:
        print("\n✓ LLM client is working correctly")
    else:
        print("\n✗ LLM client test failed - check API key and network")
