from setuptools import setup, find_packages


setup(
    name='fsm_eigenvalue',
    version='1.0.0',
    url='https://github.com/petarmaric/fsm_eigenvalue',
    license='BSD',
    author='Petar Maric',
    author_email='petarmaric@uns.ac.rs',
    description='Console app and Python API implementing a generalization of '\
                'eigenvalue problem within the harmonic coupled finite strip '\
                'method, used for parametric modeling of static and dynamic '\
                'inelastic buckling, free vibration, damage and failure in '\
                'prismatic shell structures.',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
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
