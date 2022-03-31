from setuptools import setup, find_packages


def main():

    setup(name='cell_cycle_gating',
          version='0.1.0',
          description='Deep Dye Drop based cell cycle gating',
          author='Marc Hafner, Kartik Subramanian',
          author_email='marc_hafner@hms.harvard.edu',
          url='http://github.com/datarail/DrugResponse',
          packages=find_packages(),
          install_requires=['numpy', 'pandas', 'peakutils', 'seaborn', 'matplotlib'],
          keywords=['systems', 'biology', 'data', 'array', 'matrix'],
          classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Visualization',
            ],
          )


if __name__ == '__main__':
    main()
