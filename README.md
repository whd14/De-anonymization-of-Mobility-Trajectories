# De-anonymization-of-Mobility-Trajectories

This software is an implementation of the following paper:

Huandong Wang, Chen Gao, Yong Li, Gang Wang, Depeng Jin, and Jingbo Sun
"De-anonymization of Mobility Trajectories: Dissecting the Gaps between Theory and Practice",
in Proceedings of the Network and Distributed System Security Symposium (NDSS), 2018.

The goal of this software is to de-aonymize ISP trajectories using external information, e.g., users' check-ins at online social network. In this software, our proposed algorithms and 9 baselines algorithms are all implemented to de-anonymize ISP trajectories.


Algorithms
----------
* GM [NDSS 2018]
* POIS [WWW 2016]
* ME [AIHC 2016]
* HMM [IEEE SP 2011]
* WYCI [WOSN 2014]
* HIST [TIFS 2016]
* NFLX [IEEE SP 2008]
* MSQ [TON 2013]
* LRCF [WWW 2013]
* MKV [WPES 2008]


Parameters
----------
Parameter    | Definition|Default|Remarks|
:------------:|:---------:|:---------:|:---------:|
|input_path|Path of trace data|./input|-|
|output_path|Path of output|./output|-|
|save_weights|Save candidates' weight or not|False|Weight will be saved in output_path if set to True|
|region_centre|Location of regions' centers|./RegionCenters|RegionID and its location|
|region_popularity|Popularity of regions|./RegionPopularity|Only used in **LRCF** & **POIS**|
|Tmax|The number of time-bins|24\*8|Only used in **HMM** & **GMM**|
|scale|Spatial granularity|0.01|Granularity of lng&lat|
|K|Top-K in hit-precision|10|-|

How to use (taking POIS as an example)
----------
1. Get help
```bash
$ python algorithms/pois.py --help
```
2. Run
```bash
$ python algorithms/pois.py \
    --input_path ./algorithms/input/ \
    --output_path ./algorithms/output/pois \
    --save_weights False \
    --region_centre ./algorithms/RegionCenters \
    --region_popularity ./algorithms/RegionPopularity \
    --scale 0.01 \
    --K 10
```

Example
----------
We leave a synthesized dataset of trajectories [data_example](https://github.com/whd14/De-anonymization-of-Mobility-Trajectories/blob/master/algorithms/input/), which can be used as input of the de-anonymization algorithms. 

In this example, the first line [data_example](https://github.com/whd14/De-anonymization-of-Mobility-Trajectories/blob/master/algorithms/input/) represents the external trajectories.

```bash
user_000001;|163,31.7324327_120.9655657|181,31.7324327_120.9655657|
```

The string before semi-colon is the userID. Spatio-temporal points in the trajectories are divided by vertical bar. The element before comma is the ID of the time-bin (1-hour in this exmaple), and the element after comma is latitude and longitude.

From the second line to the last line are ISP trajectory, of which an example is:

```bash
user_000010;|163,fe4c21d14228bd83d838ee6562332c0b|
```
Different with external trajectory, locations in ISP trajectories are represented by the ID of region. The geographic coordinates of the centers of regions can be found in [algorithm\RegionCenters](https://github.com/whd14/De-anonymization-of-Mobility-Trajectories/blob/master/algorithms/RegionCenters).

**Run the code**

You can simply use the following command to run all the algorithms on [data_example](https://github.com/whd14/De-anonymization-of-Mobility-Trajectories/blob/master/algorithms/input/):
```bash
./run.sh
```
`algorithms/output/XX/result` is the obtained result, where `XX` represents the name of the de-anonymization algorithm. For each input file,
there is a line in `result` with two fields separated by "\t". The first field is the name of input file. The second field is the hit-precision.

If you set `save_weights` as `True`, weight will be saved in output_path with the same name of input file.

**Other**

1. [RegionPopularity](https://github.com/whd14/De-anonymization-of-Mobility-Trajectories/blob/master/algorithms/RegionPopularity) gives the total number of access logs in each region. De-anonymization algorithms of LRCF and POIS need this file as input.

2. You can use the scripts `ParaLearn.py` to learn the parameters of GM algorithm.

REFERENCES
==========

Please cite the reference appropriately in case they are used.

[NDSS 2018]  H. Wang, C. Gao, Y. Li, G. Wang, D. Jin, J. Sun, "De-anonymization of Mobility Trajectories: Dissecting the Gaps between Theory and Practice," in Proc. NDSS, 2018.

[WWW 2016] C. Riederer, Y. Kim, A. Chaintreau, N. Korula, and S. Lattanzi, “Linking users across domains with location data: Theory and validation,” in Proc. WWW, 2016.

[AIHC 2016] A. Cecaj, M. Mamei, and F. Zambonelli, “Re-identification and information fusion between anonymized cdr and social network data,” Journal of Ambient Intelligence and Humanized Computing, vol. 7, no. 1, pp. 83–96, 2016.

[WOSN 2014] L. Rossi and M. Musolesi, “It’s the way you check-in: identifying users in location-based social networks,” in Proc. ACM WOSN, 2014.

[TIFS 2016] F. M. Naini, J. Unnikrishnan, P. Thiran, and M. Vetterli, “Where you are is who you are: User identification by matching statistics,” IEEE Transactions on Information Forensics and Security (TIFS), vol. 11, no. 2, pp. 358–372, 2016.

[IEEE SP 2008] A. Narayanan and V. Shmatikov, “Robust de-anonymization of large sparse datasets,” in Proc. IEEE SP, 2008.

[IEEE SP 2011] R. Shokri, G. Theodorakopoulos, J.-Y. Le Boudec, and J.-P. Hubaux, “Quantifying location privacy,” in Proc. IEEE SP, 2011.

[TON 2013] C. Y. Ma, D. K. Yau, N. K. Yip, and N. S. Rao, “Privacy vulnerability of published anonymous mobility traces,” IEEE/ACM Transactions on Networking (TON), vol. 21, no. 3, pp. 720–733, 2013.

[WWW 2013] O. Goga, H. Lei, S. H. K. Parthasarathi, G. Friedland, R. Sommer, and
R. Teixeira, “Exploiting innocuous activity for correlating users across
sites,” in Proc. WWW, 2013.

[WPES 2008] Y. De Mulder, G. Danezis, L. Batina, and B. Preneel, “Identification via
location-profiling in gsm networks,” in Proc. ACM WPES, 2008.
