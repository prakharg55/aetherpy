# aetherpy
Python package that supports Aether analysis

To use:

git clone https://github.com/AetherModel/aetherpy

cd aetherpy

git checkout develop

sudo python setup.py install

cp aetherpy/run_plot_model_results.py [to_someplace_in_your_PATH]

then you should be able to do something like:

run_plot_model_results.py -var=3 -alt=120 3DALL*.bin

from whatever directory you have files. You can see the help by typing:

run_plot_model_results.py -h

If you want do development, then you can edit the original files in the aetherpy directories (io, plot, utils).  You can then git commit things.  But, in order to actually use the changes that you implemented in the original files, you have to go into the aetherpy base directory, and then do:

sudo python setup.py install

again.  This overwrites the files in the system directory so when you run the codes again, you will se your changes implemented.
