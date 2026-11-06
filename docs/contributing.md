# Contributing to OpenMembraneBuilder

OpenMembraneBuilder is an open-source community-driven platform focused on reproducible, extensible, and collaborative scientific workflows and infrastructure.

We welcome contributions from developers, researchers, testers, and designers.

---

# Ways to Contribute

Contributions can broadly be divided into two categories:

## 1. Testing and Validation

Testing is a critical part of the project and includes:

- bug reporting
- workflow validation
- benchmarking
- documentation improvements
- reproducibility checks
- platform compatibility testing
- usability feedback

Scientific software quality strongly depends on rigorous testing and validation from the community.

---

## 2. Module and Feature Development

We also welcome technical contributions such as:

- new modules
- workflow integrations
- UI/UX improvements
- API integrations
- visualization components
- documentation systems
- deployment tooling

The project is designed to evolve modularly, and community-developed components are highly encouraged.

---

# Academic Contribution and Authorship

OpenMembraneBuilder is developed with an open-science philosophy.

Major contributors to the project may be recognized through authorship in associated academic manuscripts and preprints, depending on the scope and significance of their contributions.

## Living Preprint Model

The project aims to maintain a continuously updated preprint accompanying major releases of the platform.

For each major release:

- contributors may be acknowledged or included as co-authors
- release-specific contributions may be documented
- major technical and scientific advancements may be summarized

This creates a transparent and evolving record of community contributions.

---

# Future Journal Submission

As the platform matures and stabilizes, we plan to submit an associated manuscript to a peer-reviewed scientific journal (e.g. methods- or infrastructure-oriented journals).

The goal is to establish OpenMembraneBuilder as:

- a robust scientific platform
- a reproducible open framework
- a sustainable community-driven infrastructure

---

# Development setup

```bash
git clone -b development https://github.com/lolicato/OpenMembraneBuilder.git
cd OpenMembraneBuilder
pip install -e ".[docs]"     # editable install + documentation tools
omb gui                      # code changes apply on the next page reload
mkdocs serve                 # documentation preview on http://127.0.0.1:8000
```

The [code reference](code/index.md) explains how the code is organised. Every function has a
plain-language explanation (docstring) in Google style (`Args:`, `Returns:`); please keep that
style for new code, so it appears in this documentation.

---

# How to Get Involved

If you would like to contribute, collaborate, or discuss ideas, please contact us at:

dev@open-builder.com

Please include:

- your background/interests
- areas you would like to contribute to