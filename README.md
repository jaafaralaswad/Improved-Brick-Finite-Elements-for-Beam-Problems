![Python Version](https://img.shields.io/badge/python-3.12-blue)
![OS](https://img.shields.io/badge/os-ubuntu%20%7C%20macos%20%7C%20windows-blue)
![License](https://img.shields.io/badge/license-MIT-green)

[![codecov](https://codecov.io/gh/jaafaralaswad/Assignment-2/branch/main/graph/badge.svg)](https://codecov.io/gh/jaafaralaswad/Assignment-2) ![GitHub Actions](https://github.com/jaafaralaswad/Assignment-2/actions/workflows/tests.yml/badge.svg)

# ME700 Assignment 2

## Table of Contents

- [Introduction](#introduction)
- [Conda Environment, Installation, and Testing](#conda-environment-installation-and-testing)
- [The Direct Stiffness Method](#the-direct-stiffness-method)
- [Tutorials](#tutorials)
- [More Information](#more-information)

## Introduction
This repository presents the work developed to fulfill the requirements of Assignment 2 for the course ME700.


## Conda environment, install, and testing

This procedure is very similar to what we did in class. First, you need to download the repository and unzip it. Then, to install the package, use:

```bash
conda create --name assignment-2-env python=3.12
```

After creating the environment (it might have already been created by you earlier), make sure to activate it, use:

```bash
conda activate assignment-2-env
```

Check that you have Python 3.12 in the environment. To do so, use:

```bash
python --version
```

Create an editable install of the assignemnt codes. Use the following line making sure you are in the correct directory:

```bash
pip install -e .
```

You must do this in the correct directory; in order to make sure, replace the dot at the end by the directory of the folder "Assignment-2-main" that you unzipped earlier: For example, on my computer, the line would appear as follows:

```bash
pip install -e /Users/jaafaralaswad/Downloads/Assignment-2-main
```

Now, you can test the code, make sure you are in the tests directory. You can know in which directory you are using:

```bash
pwd
```

Navigate to the tests folder using the command:

```bash
cd
```

On my computer, to be in the tests folder, I would use:

```bash
cd /Users/jaafaralaswad/Downloads/Assignment-2-main/tests
```


Once you are in the tests directory, use the following to run the tests:

```bash
pytest -s test_main.py
```

Code coverage should be 100%.

To run the tutorial, make sure you are in the tutorials directory. You can navigate their as you navigated to the tests folder. On my computer, I would use:

```bash
cd /Users/jaafaralaswad/Downloads/Assignment-2-main/tutorials
```

Once you are there, you can use:

```bash
pip install jupyter
```

Depending on which tutorial you want to use, you should run one of the following lines:


```bash
jupyter notebook direct_stiffness_method.ipynb
```

```bash
jupyter notebook technical_correctness_check.ipynb
```

A Jupyter notebook will pop up.



## The Direct Stiffness Method
To be written


## Tutorials

To be written.

## More information

To be written
