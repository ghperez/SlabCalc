from distutils.core import setup

setup(
    name="SlabCalc",
    version="0.2",
    author="Gabriel H. Perez",
    author_email="gabriel.perez@aluno.ufabc.edu.br",
    description="",
    url="https://github.com/ghperez/slab_builder",
    packages=['slabCalc'],
    requires=["molSimplify","qe"]
    )