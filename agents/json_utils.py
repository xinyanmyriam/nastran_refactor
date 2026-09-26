"""
Shared JSON parsing utilities for all agents.

Handles common LLM output issues:
- LaTeX escape sequences (\\frac, \\(, \\))
- Truncated JSON
- Markdown code blocks
- Invalid escape characters
"""

import json
import re


def fix_invalid_escapes(s: str) -> str:
    """Remove invalid JSON escape sequences (e.g., LaTeX \\frac, \\(, \\))."""
    # Valid JSON escapes: \" \\ \/ \b \f \n \r \t \uXXXX
    return re.sub(r'\\([^"\\/bfnrtu])', r'\1', s)


def extract_json_block(response: str) -> str:
    """Extract JSON from a response that may contain markdown code blocks."""
    if "```json" in response:
        return response.split("```json")[1].split("```")[0].strip()
    elif "```" in response:
        # Check if the first code block looks like JSON
        block = response.split("```")[1].split("```")[0].strip()
        if block.startswith("{") or block.startswith("["):
            return block
    # No code block, try the whole response
    # Find the first { or [ and last } or ]
    start_obj = response.find("{")
    start_arr = response.find("[")
    if start_obj == -1 and start_arr == -1:
        return response.strip()
    if start_arr == -1 or (start_obj != -1 and start_obj < start_arr):
        # Object
        end = response.rfind("}")
        if end > start_obj:
            return response[start_obj:end+1]
    else:
        # Array
        end = response.rfind("]")
        if end > start_arr:
            return response[start_arr:end+1]
    return response.strip()


def robust_json_parse(response: str, default: dict = None) -> dict:
    """
    Parse JSON from LLM response with multiple fallback strategies.
    
    Strategies:
    1. Direct parse
    2. Extract JSON block from markdown, then parse
    3. Fix invalid escapes, then parse
    4. Truncate at last valid closing brace, then parse
    5. Regex extraction of key fields
    
    Returns parsed dict or default if all strategies fail.
    """
    if default is None:
        default = {}
    
    json_str = extract_json_block(response)
    
    # Strategy 1: Direct parse
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        pass
    
    # Strategy 2: Fix invalid escapes
    fixed = fix_invalid_escapes(json_str)
    try:
        return json.loads(fixed)
    except json.JSONDecodeError:
        pass
    
    # Strategy 3: Try to fix truncated JSON by finding last valid close
    for attempt in [json_str, fixed]:
        # Try progressively shorter substrings ending with }
        for end_pos in range(len(attempt) - 1, max(0, len(attempt) - 200), -1):
            if attempt[end_pos] == '}':
                try:
                    return json.loads(attempt[:end_pos+1])
                except json.JSONDecodeError:
                    continue
    
    # Strategy 4: Regex extraction of key-value pairs
    result = dict(default)
    # Extract string fields
    for match in re.finditer(r'"(\w+)"\s*:\s*"((?:[^"\\]|\\.)*)"', json_str):
        key, value = match.group(1), match.group(2)
        if key in result or not result:
            result[key] = value
    
    # Extract number fields
    for match in re.finditer(r'"(\w+)"\s*:\s*(\d+(?:\.\d+)?)', json_str):
        key, value = match.group(1), match.group(2)
        if '.' in value:
            result[key] = float(value)
        else:
            result[key] = int(value)
    
    # Extract array fields (simple string arrays)
    for match in re.finditer(r'"(\w+)"\s*:\s*\[(.*?)\]', json_str, re.DOTALL):
        key = match.group(1)
        items = re.findall(r'"((?:[^"\\]|\\.)*)"', match.group(2))
        if items:
            result[key] = items
    
    # Extract boolean fields
    for match in re.finditer(r'"(\w+)"\s*:\s*(true|false)', json_str, re.IGNORECASE):
        result[match.group(1)] = match.group(2).lower() == 'true'
    
    return result if result != default else default
