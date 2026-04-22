"""LLM adapter for reporter query enhancement."""

from __future__ import annotations

import json
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


@dataclass
class LLMResponse:
    """LLM response structure."""

    answer: str
    confidence: str
    reasoning: list[str]
    raw_response: dict[str, Any] | None = None


class LLMAdapter(ABC):
    """Abstract base class for LLM adapters."""

    @abstractmethod
    def generate_answer(
        self,
        question: str,
        context: dict[str, Any],
        query_type: str,
    ) -> LLMResponse:
        """Generate answer using LLM."""
        pass


class OpenAIAdapter(LLMAdapter):
    """OpenAI API adapter."""

    def __init__(self, api_key: str | None = None, model: str = "gpt-4o-mini"):
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        self.model = model
        if not self.api_key:
            raise ValueError("OpenAI API key not provided")

    def generate_answer(
        self,
        question: str,
        context: dict[str, Any],
        query_type: str,
    ) -> LLMResponse:
        """Generate answer using OpenAI API."""
        try:
            import openai
        except ImportError:
            raise ImportError("openai package not installed. Run: pip install openai")

        client = openai.OpenAI(api_key=self.api_key)

        system_prompt = self._build_system_prompt(query_type)
        user_prompt = self._build_user_prompt(question, context)

        response = client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.3,
            max_tokens=800,
        )

        content = response.choices[0].message.content
        return self._parse_response(content, response.model_dump())

    def _build_system_prompt(self, query_type: str) -> str:
        """Build system prompt based on query type."""
        base = """You are Immunex Assistant, a scientific AI specialized in TCR-pMHC molecular dynamics analysis.

Your role:
- Provide concise, accurate answers based on analysis data
- Use scientific terminology appropriately
- Cite evidence from the provided data
- Be honest about limitations and uncertainties

Output format:
- Start with a direct answer (1-2 sentences)
- Follow with 2-4 evidence bullets
- End with confidence level (high/medium/low)
"""

        type_specific = {
            "cluster_status": "Focus on interface clustering states and dominant conformations.",
            "quality_status": "Focus on trajectory stability, RMSD convergence, and simulation quality.",
            "interface_summary": "Focus on buried surface area, contact patterns, and interface characteristics.",
            "hotspot_summary": "Focus on RRCS scores, key residue pairs, and interaction hotspots.",
            "design_hint": "Focus on mutation candidate sites based on hotspot analysis.",
            "top_pair": "Focus on the strongest RRCS hotspot pairs.",
            "top_residue": "Focus on the most critical residues by RRCS score.",
        }

        return base + "\n" + type_specific.get(query_type, "")

    def _build_user_prompt(self, question: str, context: dict[str, Any]) -> str:
        """Build user prompt with question and context."""
        context_str = json.dumps(context, indent=2, ensure_ascii=False)
        return f"""Question: {question}

Available data:
{context_str}

Please provide:
1. A direct answer (1-2 sentences)
2. Evidence bullets (2-4 items)
3. Confidence level (high/medium/low)

Format your response as:
ANSWER: [your answer]
EVIDENCE:
- [evidence 1]
- [evidence 2]
CONFIDENCE: [high/medium/low]
"""

    def _parse_response(self, content: str, raw: dict) -> LLMResponse:
        """Parse LLM response into structured format."""
        lines = content.strip().split("\n")
        answer = ""
        evidence = []
        confidence = "medium"

        current_section = None
        for line in lines:
            line = line.strip()
            if line.startswith("ANSWER:"):
                current_section = "answer"
                answer = line.replace("ANSWER:", "").strip()
            elif line.startswith("EVIDENCE:"):
                current_section = "evidence"
            elif line.startswith("CONFIDENCE:"):
                confidence = line.replace("CONFIDENCE:", "").strip().lower()
                current_section = None
            elif line.startswith("-") and current_section == "evidence":
                evidence.append(line.lstrip("- ").strip())
            elif current_section == "answer" and line:
                answer += " " + line

        return LLMResponse(
            answer=answer or "Unable to generate answer from available data.",
            confidence=confidence,
            reasoning=evidence,
            raw_response=raw,
        )


class AnthropicAdapter(LLMAdapter):
    """Anthropic Claude API adapter."""

    def __init__(self, api_key: str | None = None, model: str = "claude-3-5-sonnet-20241022"):
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        self.model = model
        if not self.api_key:
            raise ValueError("Anthropic API key not provided")

    def generate_answer(
        self,
        question: str,
        context: dict[str, Any],
        query_type: str,
    ) -> LLMResponse:
        """Generate answer using Anthropic API."""
        try:
            import anthropic
        except ImportError:
            raise ImportError("anthropic package not installed. Run: pip install anthropic")

        client = anthropic.Anthropic(api_key=self.api_key)

        system_prompt = self._build_system_prompt(query_type)
        user_prompt = self._build_user_prompt(question, context)

        message = client.messages.create(
            model=self.model,
            max_tokens=800,
            temperature=0.3,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}],
        )

        content = message.content[0].text
        return self._parse_response(content, message.model_dump())

    def _build_system_prompt(self, query_type: str) -> str:
        """Build system prompt (same as OpenAI)."""
        return OpenAIAdapter(api_key="dummy")._build_system_prompt(query_type)

    def _build_user_prompt(self, question: str, context: dict[str, Any]) -> str:
        """Build user prompt (same as OpenAI)."""
        return OpenAIAdapter(api_key="dummy")._build_user_prompt(question, context)

    def _parse_response(self, content: str, raw: dict) -> LLMResponse:
        """Parse response (same as OpenAI)."""
        return OpenAIAdapter(api_key="dummy")._parse_response(content, raw)


def create_llm_adapter(provider: str = "openai", **kwargs) -> LLMAdapter:
    """Factory function to create LLM adapter."""
    adapters = {
        "openai": OpenAIAdapter,
        "anthropic": AnthropicAdapter,
    }

    if provider not in adapters:
        raise ValueError(f"Unknown provider: {provider}. Available: {list(adapters.keys())}")

    return adapters[provider](**kwargs)
