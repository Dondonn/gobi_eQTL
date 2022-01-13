# Homework 02

## Submission
Zip the entire directory (the root directory in which `src`, `tests`, `environment.yml` and `README.md`) are located.
Name the archive according to the following scheme:
```
f'ex_02_{student.first_name}_{student.last_name}_{student.immatriculation_number}.zip'
```

Send the archive to `advist@rostlab.org`. Put the following in the email topic:
```
f'Exercise 02 {student.first_name} {student.last_name} {student.immatriculation_number}'
```

Submit your code before 12.11.2021 09:30.

## Environment
Create and activate the environment using the following commands:
```bash
conda env create -f environment.yml
conda activate advist-02
```

## Transcription
Complete the `transcribe` function.
Pay attention to the docstring, it contains the specification.

## Translation
Complete the `translate` function.
The codon table is provided in the module for your convenience.
Pay attention to the docstring, it contains the specification.

## Transcription and/or Translation
Complete the `nuc_to_aa` function.
Feel free to do it any way you want, but, hopefully, you will use the `transcribe` and `translate` functions.

## NGS9000
Complete the `NGS9000` class.

## Testing your code
The following code runs all the tests:

```bash
python -m pytest
```

The tests have been split into several files for your convenience. Run specific tests using the following code:
```bash
python -m pytest tests/test_transcribe.py
python -m pytest tests/test_translate.py
python -m pytest tests/test_nuc_to_aa.py
python -m pytest tests/test_ngs_9000.py
```

## Plot
Complete the `plot` function in the `src/plot.py` file. The function should generate a .png file in the `src` directory.
"# gobi_eQTL" 
