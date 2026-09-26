"""
Multi-Agent Framework for NASTRAN Legacy Code Refactoring.

This package implements a multi-agent system that automates the modernization
of NASTRAN-95 Fortran 77 code into modern Python, with equivalence verification.

Agents:
    - ComprehensionAgent: Understands Fortran source code intent
    - ArchitectureAgent: Recovers module dependencies and data flow
    - PhysicsSpecAgent: Extracts structured physics specifications
    - ModernDesignAgent: Designs modern Python architecture
    - CodeGenAgent: Generates Python implementation
    - EquivalenceAgent: Verifies behavioral equivalence
    - RefactoringOrchestrator: Coordinates the full pipeline
"""

from .orchestrator import RefactoringOrchestrator

__all__ = ["RefactoringOrchestrator"]
