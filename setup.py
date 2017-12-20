from setuptools import setup, find_packages


setup(
    name='fsm_eigenvalue',
    version='0.1.0.dev',
    url='https://bitbucket.org/petar/fsm_eigenvalue',
    license='BSD',
    author='Petar Maric',
    author_email='petarmaric@uns.ac.rs',
    description='Console app and Python API implementing a generalization of '\
                'eigenvalue problem within the harmonic coupled finite strip '\
                'method, used for parametric modeling of static and dynamic '\
                'inelastic buckling, free vibration, damage and failure in '\
                'prismatic shell structures.',
    long_description=open('README').read(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    platforms='any',
    packages=find_packages(),
    entry_points={
        'console_scripts': ['fsm_eigenvalue=fsm_eigenvalue.shell:main'],
    },
    install_requires=open('requirements.txt').read().splitlines(),
)
