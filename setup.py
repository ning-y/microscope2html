from setuptools import setup, find_packages

setup(
    name="microscope2html",
    version="1.1.0",
    description="Stitch microscope tiles and generate a standalone HTML viewer",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "microscope2html=microscope2html.cli:main",
        ],
    },
    install_requires=[
        "numpy",
        "tifffile",
        "Pillow",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
