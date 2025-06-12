# FlowTools
FlowTools is a collection of scripts to help with flow cytometry related tasks and make your lab life easier!

## FlowJo data extractor
The FlowJo data extractor tool is a script based on the [FlowKit](https://github.com/whitews/FlowKit) library that allows you to extract a tidy data dataset directly from your FlowJo file.

### Setup
#### 1. Python setup and packages
- First, you need Python on your machine to run this code. I recommend [downloading anaconda](https://www.anaconda.com/download) for your python installation.
- Now to install the required packages
  - Open your terminal:
    - run "pip install flowkit"
    - run "conda install conda-forge::pysimplegui"
#### 2. Running the script
- Open your terminal
- run "python /path/to/script.py" (on Mac you can drag and drop the script file after writing "python " in the terminal and it will add the right path for you) \
A menu will pop up asking you for the required inputs.


### Some things to note
- FlowJo bins data points before computing statistics, so the exact numbers will be slighlty different from this script. Here every single event is counted and the statistics are computed on the raw data. Arguably more accurate but should not make a difference in most cases.
- The geometric mean is computed as follows: \
```math
\left(\prod _{i=1}^{n}x_{i}\right)^{\frac {1}{n}}={\sqrt[{n}]{x_{1}x_{2}\cdots x_{n}}}
```
this only makes sense if $`x_{1}x_{2}\cdots x_{n}`$ are all greater than $0$. Due to compensation this is not the case so to make sure a geometric mean can be computed the data is transformed as such: \
$$\bar{x_i}=x_i+\min(x_{1}x_{2}\cdots x_{n})+1$$ \
the final effect is to translate all datapoints to be positive and add one to avoid having $0$ in the dataset