setup(
    name="mrsoftware",
    version=VERSION,
    description='CSE185 Final Project',
    author="Johnathan Narita and Tasha Nguyen",
    author_email="tpn006@ucsd.edu",
    packages=find_package(),
    entry_points={
        "console_scripts": [
            "mrsoftware=mrsoftware.mrsoftware:main"
        ],
    },
)