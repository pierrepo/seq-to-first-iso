default: help


env: ## Create conda env
	conda env create -f environment.yml
.PHONY: env

test: ## Run tests
	py.test tests
.PHONY: test

test-coverage: ## Run tests with coverage
	py.test --cov --cov-config .coveragerc
.PHONY: test-coverage

lint: ## Lint code
	pycodestyle seq_to_first_iso \
	&& pydocstyle seq_to_first_iso \
	&& pylint seq_to_first_iso
.PHONY: lint

compile: ## Compile for PyPI
	python setup.py bdist_wheel
.PHONY: compile

upload-to-pypi: ## Upload to PyPI
	twine upload dist/*
	# clean compiled
	rm -f dist/*.tar.gz dist/*.whl dist/*.egg
.PHONY: upload-to-pypi

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
