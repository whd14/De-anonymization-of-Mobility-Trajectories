Learning the parameters of GM algorithm
----------

**Format of input data**

Input file: `Mismatches` 

Each line represents the distance between a record in an external trajectory and
its time-adjacent records in ISP trajectory belonging to the same user. The length of time window is set to be 24 hour in this example. Thus, each line contains 24 numbers. We give 6349 records (lines) in the file `Mismatches` as an example. 

**Run the code**

You can use the following command to run the algorithms on the file `Mismatches`:
```bash
python ParaLearn.py
```
`para` is the obtained parameters. 

In addition, parameters can also be set based on empirical distribution of spatio-temporal mismatches. For example, *pi* is set to be exponential or power-law distribution, while *sigma2* can be set as 0.5^2, which represent spatial mismatches of 500 meters.
