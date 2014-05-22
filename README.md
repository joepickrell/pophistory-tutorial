#Tutorial on inferring history from genetics
## Joe Pickrell

#Introduction

The goal of this tutorial is to expose you to some of the important tools for learning about history from genetics. The standard setup for someone interested in these problems is to genotype thousands or millions of SNPs in a population of interest and then ask: what is the history of this population? The tools you will run in this tutorial provide summaries of the data that are informative about this question. However, be aware that there is no ``black box" solution--all methods can be misleading, and the main difficulty is not running software, but instead interpreting the results.

There are two example datasets in this tutorial, both taken from human population genetics. The first example is designed to show how even simple analyses can have complex interpretations, and the second is designed to give you some exposure to several of the main analysis tools used in papers about population history.

# Preliminaries

We assume you are comfortable working from the command line and have installed:
*`smartpca` (in the [EIGENSOFT package](http://www.hsph.harvard.edu/alkes-price/software/))
*treemix
*STRUCTURE
 
# Example 1
There are two directories in this tutorial, `example1/` and `example2/`. We're going to start with example 1
```
cd example1/

ls
```


\noindent In this directory is a single file in TreeMix format, \texttt{example1.treemix.gz}. Look at the first few lines of this file:
\\
\\
\texttt{zcat example1.treemix.gz $|$ head}
\\
\\
\noindent The first line contains the names of five human populations: French, Sardinian, Karitiana (a Native American population), Yoruba (Nigeria), and Han (China). The remainder of the lines contain SNP data: each entry is of the form [count1],[count2], where the counts are the numbers of observed reference and alternative alleles in a set of genotyped individuals from the population. A bioinformatics question:

\begin{enumerate}
\item Using the command line, count the number of SNPs included in this file. Hint: if you're not familiar with the unix command \texttt{wc}, run:
\\
\\
\texttt{man wc}
\end{enumerate}

\noindent Now run four-population tests on all five possible combinations of populations. To do this, use the \texttt{fourpop} command from the TreeMix package. First get the help page by running:
\\
\\
\texttt{fourpop}\\
\\
\noindent You should see that you need to give the program an input file as well as a setting for the numbers of SNPs to use in the block jackknife. This parameter ensures that LD structure does not lead to ``double-counting" of SNPs. In these data, a setting of \texttt{k} of about 500 should be sufficient. Now run \texttt{fourpop}:
\\
\\
\texttt{fourpop -i example1.treemix.gz -k 500}\\
\\
\noindent The program will print to the screen all possible four-population trees in the data in the format:
\\
\\
\noindent [tested tree] [$f_4$ statistic] [standard error] [Z-score]
\\
\\
\noindent In the tested tree, the format is A,B;C,D, such that the computed $f_4$ statistic is (A-B)(C-D). Recall that a Z-score with an absolute value of 3 or greater corresponds to a P-value of 0.001. Questions:

\begin{enumerate}
\item What is the best tree for the populations: Yoruba, Han, French, Karitiana?
\item Does this tree fit the data? If not, how do you interpret the result?
\end{enumerate}

\noindent Another piece of information that might be useful for interpreting these data comes from $f_3$ statistics. These can be computed using the program \texttt{threepop}:
\\
\\
\texttt{threepop -i example1.treemix.gz -k 500}
\\
\\
\noindent Now the output is of the form:
\\
\\
\noindent [tested tree][$f_3$ statistic][standard error][Z-score]
\\
\\
\noindent In the tested tree, the format is A;B,C, such that the computed $f_3$ statistic is (A-B)(A-C). Recall that for $f_3$ statistics, only significantly \emph{negative} statistics are informative. Questions:
\begin{enumerate}
\item Are there any significantly negative $f_3$ statistics in these data? 
\item How do you interpret the combination of these $f_3$ and $f_4$ statistics?
\end{enumerate}


\section{Example 2}
Now let's go to the second example:
\\
\\
\texttt{cd $\sim$/ACAD/Wed\_Pickrell/example2/}
\\
\texttt{ls}
\\
\\
\noindent In this directory is data from five human populations in three different formats. The five populations are a northern European population (population code: CEU), Chinese (CHB), Yoruba from Nigeria (YRI), Maasai from Kenya (MKK), and Gujarati from India (GIH). First is STRUCTURE format: \texttt{example2.struct\_in}. Then is EIGENSTRAT format: \texttt{example2.eigenstratgeno}, \texttt{example2.snp}, and \texttt{example2.ind}. Finally is TreeMix format: \texttt{example2.treemix.gz}. Take some time to look over the file formats. For example, look at the format for TreeMix:
\\
\\
\texttt{zcat example2.treemix.gz $|$ head}
\\
\\
\noindent The file format is simple an ordered list of the counts of the alleles of SNPs in five populations, with the population codes listed in the header. In the other file formats are the genotypes for the individuals at different sets of SNPs.
 
We now want to understand the relationships between the five populations in these data. The goal here is to see how different approaches to looking at the data tell us similar or different things about history.
\subsection{STRUCTURE}
In the file \texttt{example2.struct\_in} is the input file for the program \texttt{STRUCTURE}. It contains 75 individuals from the five populations genotyped at 3,000 SNPs. In this file, populations have numbers instead of names: CEU is population 1, CHB is population 2, GIH is population 3, MKK is population 4, and YRI is population 5. Start STRUCTURE by moving into its directory and starting the GUI:
\\
\\
\texttt{cd $\sim$/frontend/}\\
\texttt{./structure}
\\
\\
\noindent In the File menu, select ``New project" and follow the instructions. Most of the steps involve telling the program the format of the file: check the boxes to tell the program that there is a ``row of marker names" and ``map distances between loci", and that the file is in the ``special format". Then check the boxes for ``Individual ID for each individual" and ``Putative population origin for each individual". The data is not phased.

After the file has loaded, go to the ``Parameter set" menu and select ``New". Parameter estimation is done by MCMC, so the program needs to know how long to run and how many iterations to discard. For our purposes here, set the burn-in to 1,000 iterations, at a few thousand additional iterations should be enough to get an idea of how the data look. 

After creating the parameter set, return to the ``Parameter set" menu and select ``Run". The program will ask the number of populations--enter 2 to begin. After the run has completed, click on the results file and play with the different visualization methods. Questions:

\begin{enumerate}
\item Run STRUCTURE with different settings for the number of populations, how do the results change?
\item How do you interpret these results?
\end{enumerate}
\subsection{PCA}
In the files \texttt{example2.eigenstratgene}, \texttt{example2.snp}, and \texttt{example2.ind} are 116,565 SNPs genotypes from 675 individuals from the five populations. A quick bioinformatics question:

\begin{enumerate}
\item Using R, count the number of individuals per population. Hint: in R, look at the \texttt{table} function by running:
\\
\\ 
\texttt{R}\\
\texttt{>read.table("example2.ind", as.is = T)}\\
\texttt{>?table}\\
\texttt{>q() }[to return to the command line]\\
\end{enumerate}

\noindent Now return to the command line and open the file \texttt{par.example2}. This is the parameter file for \texttt{smartpca}. To tell the program to use these parameters, run:
\\
\\
\texttt{smartpca -p par.example2}
\\
\\
\noindent The program will output the positions of each individual on the major axes of variation in the sample into \texttt{example.evec}. This file allows for visualization of the population structure in the sample. The best way to visualize the output is using R:
\\
\\
\texttt{R}\\
\texttt{>d = read.table("example.evec", as.is = T)}\\
\texttt{>plot(d[,2], d[,3], xlab = "PC1", ylab = "PC2")}\\
\\
\noindent In these plots, each individual is a point. It would be nice to color each individual according to their population of origin. Luckily this information is in the seventh column of the output file. To color the Yoruba (YRI) population in blue, continue in R:
\\
\\
\texttt{>tmp = d[d[,7] == "YRI",]}\\
\texttt{>points(tmp[,2], tmp[,3], col = "blue", pch = 20)}\\
\\
\\
\noindent Questions:
\begin{enumerate}
\item Produce a PCA plot with all five populations plotted in different colors. Do this for PC1 versus PC2, as well as PC2 versus PC3. 
\item How do you interpret this plot?
\end{enumerate}

\subsection{TreeMix}
TreeMix is a program for building trees and graphs of populations. The input file is in \texttt{example2.treemix.gz}. Again examine the input file:
\\
\\
\texttt{zcat example2.treemix.gz $|$ head}
\\
\\
\noindent To run TreeMix, as for three- and four-population test, you need to set a window size. Again we'll set it at 500. It's also often useful to set an outgroup using the \texttt{-root} flag:
\\
\\
\texttt{treemix -i example2.treemix.gz -k 500 -root YRI}
\\
\\
\noindent TreeMix will now output a number of files with the stem TreeMix.*. The tree can be plotted in R:
\\
\\
\texttt{R}\\
\texttt{>source("plotting\_funcs.R")}\\
\texttt{>plot\_tree("TreeMix")}\\
\\
\\
\noindent It's also useful to have a measure of how well this tree fits the data. To plot the residual fit of the model, while still in R (and having loaded the functions in \texttt{plotting\_funcs.R}, run:
\\
\\
\texttt{>plot\_resid("TreeMix", "pop\_order")}\\
\\
\noindent Each entry in the displayed matrix shows how well the model accounts for the observed relationship between the pair of populations. Questions:
\begin{enumerate}
\item How do you interpret the tree and residual plot?
\end{enumerate}

\noindent TreeMix can also add migration to a tree. To tell TreeMix the number of migration events to add, use the \texttt{-m} flag:
\\
\\
\texttt{treemix -i example2.treemix.gz -k 500 -root YRI -m 1}
\\
\\
\noindent Questions:
\begin{enumerate}
\item Plot the graph and residuals using the same commands as before. How do you interpret them?
\item Run TreeMix with different parameter settings (of \texttt{-root}, \texttt{-k}, and \texttt{-m}. How robust are your results to these choices?
\end{enumerate}

\subsection{Three- and four-population tests}
As in the first example, you can also calculate $f_3$ and $f_4$ statistics on the file in TreeMix format. Questions:

\begin{enumerate}
\item How do you interpret the $f$-statistics?
\item Putting everything together, what is a historical scenario for these populations that is consistent with the STRUCTURE, PCA, TreeMix, and $f$-statistic results?
\end{enumerate}

