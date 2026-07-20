from setuptools import setup, find_packages

setup(
    name="mevd_grn",
    version="0.1.0",
    description="Multi-Evidence Distillation for Gene Regulatory Network Inference",
    packages=find_packages(where="."),
    python_requires=">=3.10",
)
