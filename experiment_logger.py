"""
Experiment Logger for NASTRAN Refactoring.

Provides comprehensive, timestamped logging for all experiments.
Every LLM call, agent action, quality gate decision, and error is recorded.
Failed experiments are preserved with full context for debugging.
"""

import json
import logging
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any


class ExperimentLogger:
    """
    Comprehensive experiment logger with timestamped file and console output.
    
    Creates a structured log directory for each experiment run:
        results/<experiment_id>/
        ├── experiment.log          # Full text log (human-readable)
        ├── experiment_meta.json    # Experiment metadata
        ├── llm_calls/              # Raw LLM request/response pairs
        │   ├── call_001.json
        │   ├── call_002.json
        │   └── ...
        ├── artifacts/              # Intermediate artifacts (JSON)
        │   ├── code_intents.json
        │   ├── architecture.json
        │   └── ...
        └── generated_code.py       # Final generated code (if any)
    """
    
    def __init__(self, experiment_id: str, base_dir: str = "results"):
        """
        Args:
            experiment_id: Unique identifier for this experiment run.
            base_dir: Base directory for all experiment results.
        """
        self.experiment_id = experiment_id
        self.start_time = datetime.now()
        self.start_timestamp = time.time()
        
        # Create directory structure
        self.base_dir = Path(base_dir)
        self.run_dir = self.base_dir / experiment_id
        self.run_dir.mkdir(parents=True, exist_ok=True)
        (self.run_dir / "llm_calls").mkdir(exist_ok=True)
        (self.run_dir / "artifacts").mkdir(exist_ok=True)
        
        # Setup Python logger
        self.logger = logging.getLogger(f"nastran_refactor.{experiment_id}")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        
        # File handler - full detail
        log_file = self.run_dir / "experiment.log"
        fh = logging.FileHandler(log_file, encoding='utf-8', mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter(
            '%(asctime)s [%(levelname)s] %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        ))
        self.logger.addHandler(fh)
        
        # Console handler - info and above.
        # The Windows console defaults to cp1252, which cannot encode the
        # box-drawing and check-mark characters used in stage banners. That
        # raised a UnicodeEncodeError on almost every log line, burying real
        # errors in traceback noise. Reconfigure the stream to UTF-8 and fall
        # back to replacing unencodable characters rather than raising.
        stream = sys.stdout
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            # Older streams, or a stream that does not support reconfigure:
            # wrap it so logging can never crash on encoding.
            import io
            try:
                stream = io.TextIOWrapper(
                    sys.stdout.buffer, encoding="utf-8", errors="replace",
                    line_buffering=True)
            except (AttributeError, ValueError):
                stream = sys.stdout

        ch = logging.StreamHandler(stream)
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter(
            '%(asctime)s [%(levelname)s] %(message)s',
            datefmt='%H:%M:%S'
        ))
        self.logger.addHandler(ch)
        
        # Counters
        self.llm_call_count = 0
        self.stage_timings: dict[str, dict] = {}
        self.errors: list[dict] = []
        self.quality_gates: dict[str, dict] = {}
        
        # Write initial metadata
        self._write_meta({
            "experiment_id": experiment_id,
            "start_time": self.start_time.isoformat(),
            "status": "running",
        })
        
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Experiment started: {experiment_id}")
        self.logger.info(f"Time: {self.start_time.isoformat()}")
        self.logger.info(f"Output directory: {self.run_dir}")
        self.logger.info(f"{'='*60}")
    
    def log_stage_start(self, stage_name: str):
        """Log the start of a pipeline stage."""
        now = time.time()
        self.stage_timings[stage_name] = {"start": now, "start_iso": datetime.now().isoformat()}
        self.logger.info(f"")
        self.logger.info(f"{'-'*60}")
        self.logger.info(f"STAGE: {stage_name}")
        # ASCII only. A U+2500 box-drawing character here rendered as mojibake
        # ("Gamma-o-C" runs) whenever the log was read back under cp1252.
        self.logger.info(f"{'-'*60}")
    
    def log_stage_end(self, stage_name: str, success: bool, details: str = ""):
        """Log the end of a pipeline stage."""
        now = time.time()
        if stage_name in self.stage_timings:
            elapsed = now - self.stage_timings[stage_name]["start"]
            self.stage_timings[stage_name]["end"] = now
            self.stage_timings[stage_name]["end_iso"] = datetime.now().isoformat()
            self.stage_timings[stage_name]["elapsed_seconds"] = round(elapsed, 2)
            self.stage_timings[stage_name]["success"] = success
        
        # ASCII only: log lines are parsed by the analysis scripts and read in
        # consoles whose codepage is not UTF-8.
        status = "SUCCESS" if success else "FAILED"
        self.logger.info(f"  Stage [{stage_name}] completed: {status} ({elapsed:.1f}s)")
        if details:
            self.logger.info(f"  Details: {details}")
    
    def log_quality_gate(self, stage_name: str, passed: bool, reason: str = ""):
        """Log a quality gate decision."""
        self.quality_gates[stage_name] = {
            "passed": passed,
            "reason": reason,
            "timestamp": datetime.now().isoformat(),
        }
        status = "PASS" if passed else "FAIL"
        self.logger.info(f"  Quality Gate [{stage_name}]: {status}")
        if reason:
            self.logger.info(f"    Reason: {reason}")
    
    def log_llm_call(
        self,
        agent_name: str,
        system_prompt: str,
        user_prompt: str,
        response: str,
        elapsed_seconds: float,
        tokens_used: int = 0,
        error: str = None,
    ):
        """Log a single LLM API call with full request/response."""
        self.llm_call_count += 1
        call_id = f"call_{self.llm_call_count:03d}"
        
        call_record = {
            "call_id": call_id,
            "timestamp": datetime.now().isoformat(),
            "agent": agent_name,
            "elapsed_seconds": round(elapsed_seconds, 2),
            "tokens_used": tokens_used,
            "system_prompt_length": len(system_prompt),
            "user_prompt_length": len(user_prompt),
            "response_length": len(response) if response else 0,
            "error": error,
            "system_prompt": system_prompt,
            "user_prompt": user_prompt,
            "response": response,
        }
        
        # Save to file
        call_file = self.run_dir / "llm_calls" / f"{call_id}.json"
        call_file.write_text(json.dumps(call_record, indent=2, ensure_ascii=False), encoding='utf-8')
        
        # Log summary
        if error:
            self.logger.warning(
                f"  LLM [{call_id}] {agent_name}: ERROR after {elapsed_seconds:.1f}s - {error}"
            )
        else:
            self.logger.debug(
                f"  LLM [{call_id}] {agent_name}: {len(response)} chars, "
                f"{elapsed_seconds:.1f}s, {tokens_used} tokens"
            )
    
    def log_artifact(self, name: str, data: Any, is_json: bool = True):
        """Save an intermediate artifact atomically."""
        if is_json:
            filepath = self.run_dir / "artifacts" / f"{name}.json"
            content = json.dumps(
                data, indent=2, default=str, ensure_ascii=False)
        else:
            filepath = self.run_dir / "artifacts" / name
            content = str(data)

        # A process interruption must leave either the previous complete file
        # or the new complete file, never a partially written JSON record.
        temp_path = filepath.with_name(filepath.name + ".tmp")
        temp_path.write_text(content, encoding="utf-8")
        temp_path.replace(filepath)

        self.logger.debug(f"  Artifact saved: {filepath.name}")
    
    def log_error(self, stage: str, error: Exception, context: str = ""):
        """Log an error with full traceback."""
        error_record = {
            "timestamp": datetime.now().isoformat(),
            "stage": stage,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "traceback": traceback.format_exc(),
            "context": context,
        }
        self.errors.append(error_record)
        
        self.logger.error(f"  ERROR in [{stage}]: {type(error).__name__}: {error}")
        self.logger.debug(f"  Traceback:\n{traceback.format_exc()}")
        
        # Also save error to file
        error_file = self.run_dir / "artifacts" / f"error_{stage}_{len(self.errors)}.json"
        error_file.write_text(json.dumps(error_record, indent=2, ensure_ascii=False), encoding='utf-8')
    
    def log_info(self, message: str):
        """Log an informational message."""
        self.logger.info(f"  {message}")
    
    def log_debug(self, message: str):
        """Log a debug message (file only)."""
        self.logger.debug(f"  {message}")
    
    def finalize(self, success: bool, summary: str = ""):
        """Finalize the experiment log."""
        end_time = datetime.now()
        total_elapsed = time.time() - self.start_timestamp
        
        self.logger.info(f"")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"EXPERIMENT COMPLETE")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"  Status: {'SUCCESS' if success else 'FAILED'}")
        self.logger.info(f"  Total time: {total_elapsed:.1f}s")
        self.logger.info(f"  LLM calls: {self.llm_call_count}")
        self.logger.info(f"  Errors: {len(self.errors)}")
        self.logger.info(f"  Quality gates: {self.quality_gates}")
        if summary:
            self.logger.info(f"  Summary: {summary}")
        self.logger.info(f"  Results: {self.run_dir}")
        self.logger.info(f"{'='*60}")
        
        # Update metadata
        self._write_meta({
            "experiment_id": self.experiment_id,
            "start_time": self.start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "total_elapsed_seconds": round(total_elapsed, 2),
            "status": "success" if success else "failed",
            "summary": summary,
            "llm_calls": self.llm_call_count,
            "errors_count": len(self.errors),
            "stage_timings": self.stage_timings,
            "quality_gates": self.quality_gates,
            "errors": self.errors,
        })
    
    def _write_meta(self, meta: dict):
        """Atomically write experiment metadata."""
        meta_file = self.run_dir / "experiment_meta.json"
        temp_file = meta_file.with_name(meta_file.name + ".tmp")
        temp_file.write_text(
            json.dumps(meta, indent=2, default=str, ensure_ascii=False),
            encoding="utf-8"
        )
        temp_file.replace(meta_file)


def create_experiment_id(model_name: str, phase: str = "phase1") -> str:
    """Create a unique experiment ID with timestamp."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{phase}_{model_name}_{timestamp}"
