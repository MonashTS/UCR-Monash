# UCR-Monash
Monash UCR suite implementation showcasing EAPrunedDTW.
This page support our paper **Early Abandoning PrunedDTW and its application to similarity search**
available on arxiv https://arxiv.org/abs/2010.05371 (submitted to SDM2021).

We developed EAPruneDTW, a version of PrunedDTW suited to early abandon.
You can find the source code of PruneDTW on its website http://sites.labic.icmc.usp.br/prunedDTW/ .

PrunedDTW was later use in the UCR suite ( **Searching and Mining Trillions ofTime SeriesSubsequencesunder Dynamic Time Warping** https://www.cs.ucr.edu/~eamonn/UCRsuite.html ),
leading to the UCR USP suite (**Speeding Up Similarity Search Under Dynamic Time Warping by Pruning Unpromising Alignments** https://sites.google.com/view/ucruspsuite ).

We implemented EAPrunedDTW in the UCR suite, leading to our own UCR MON suite.
To be clear, the UCR USP and UCR MON suites are directly derived from UCR suite: the code is the same except for the implementation of the dtw function.
The UCR suite uses lower bound to reduce computation time.
We also create a UCR MON nolb which only uses EAPrunedDTW an no lower bound at all.

We ran the same experiments as in the UCR USP paper on AMD opteron 6338P at 2.3Ghz compiled with g++ 9.3.0, with -O3 flag.
There are 6 datasets, and for each datasets we run an experiment over the combination of  5 queries,
5 different window ratios (0.1, 0.2, 0.3, 0.4, 0.5) and 4 lengths (128, 256, 512, 1024) for a total of 100 experiments per dataset.
Hence, we have 600 results per implementations, 2400 in total (for UCR, UCR USP, UCR MON, UCR MON nolb).

This repository contains:
  *  Our executables
  *  Our results
  *  Our shell script to download the datasets from the UCR USP website
  *  Our python script to generate the plots.


# How-to
The following instructions should work on most Unix/Linux platforms.

Start by downloading the datasets by launching the script in a terminal.
It may require to be set as executable, and definitely requires `curl` and `unzip`.
If you have the `parallel` tool, it will unzip the dataset in parallel.
```
chmod u+x dl_dataset.sh
./dl_dataset.sh
```

This will create a DATASETS directory, download and unzip the datasets.
You can inspect the script to see how it works.
It's weird because the datasets are stored on a google drive,
different download techniques are required depending on their size... yup!
Last but not least, the PPG datasets link is not a zip, so every queries and data are individually downloaded.

To compile a executable, navigate to its directory and launch make,
or compile it with you favourite compiler.
```
cd execs/UCR_MON
make
```


To launch a computation from the top directory (containing this README.md file),
do something like:
```
./execs/UCR_MON/UCR_MON DATASETS/FOG_LEG/Data.txt DATASETS/FOG_LEG/Query1.txt 128 0.5
```
The results will be printed on your console. You'll need a bit of redirection magic to save the result, e.g.
```
./execs/UCR_MON/UCR_MON DATASETS/FOG_LEG/Data.txt DATASETS/FOG_LEG/Query1.txt 128 0.5 > output
```
You'll still see the progress indicator (some little dot '.') as it is printed on stderr.

# About our results

The results are in the eponymous directory, organized by executable.
The python script living their can be launched to load the data.
While loading, it checks for discrepancies in the reported location and distance of the closest subsequence.
It will then produce the plot as presented in our paper, and compute:
  * the average speedup of UCR MON compare to UCR and UCR USP for all queries and window with a length of 1024,
  * check where UCR MON is slower than UCR and UCR USP

Note: we are using python 3.
If your python installation points to a python 2, you may need to use `python3 main.py` instead.

```
cd results
python main.py
```
