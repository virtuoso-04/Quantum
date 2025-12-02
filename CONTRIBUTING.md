# Quantum Drug Discovery Platform - Contributing Guide

Thank you for your interest in contributing to the Quantum Drug Discovery Platform! 🎉

## How to Contribute

### Reporting Issues
- Use GitHub Issues to report bugs or request features
- Provide detailed descriptions and reproducible examples
- Include your environment details (Python version, OS, etc.)

### Pull Requests
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add/update tests if applicable
5. Ensure code follows style guidelines (run `black` and `flake8`)
6. Commit with clear messages (`git commit -m 'Add amazing feature'`)
7. Push to your branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Code Style
- Follow PEP 8 guidelines
- Use `black` for code formatting: `black src/ examples/`
- Use type hints where possible
- Add docstrings to all public functions/classes
- Keep functions focused and modular

### Testing
- Write unit tests for new functionality
- Run existing tests: `pytest tests/`
- Aim for high code coverage

### Documentation
- Update README.md if adding new features
- Add docstrings with examples
- Update tutorials if changing API

## Priority Areas

### High Priority
- [ ] Comprehensive test suite (unit tests, integration tests)
- [ ] Enhanced solvation models (GBSA, PBSA)
- [ ] Molecular dynamics integration
- [ ] Real quantum hardware support (IBM Quantum, AWS Braket)

### Medium Priority
- [ ] Larger basis set support (cc-pVDZ, cc-pVTZ)
- [ ] Additional ansätze (ADAPT-VQE, QEB)
- [ ] Machine learning integration (energy prediction)
- [ ] PDB file parsing for protein structures
- [ ] Automated geometry optimization

### Nice to Have
- [ ] Web dashboard with real-time updates
- [ ] Cloud deployment options
- [ ] More visualization options (3D molecular viewer)
- [ ] Benchmarking suite
- [ ] CI/CD pipeline

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/Quantum.git
cd Quantum

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install development dependencies
pip install -r requirements.txt
pip install pytest black flake8 mypy

# Run tests
pytest tests/

# Format code
black src/ examples/
flake8 src/ examples/
```

## Questions?

Feel free to open an issue or reach out via email!

Happy coding! 🚀
