# OPUS-Fold

<img src="./images/figure1.png" width="80%" height="80%"/>

We propose a protein folding framework, named OPUS-Fold, which can integrate various methods for subproblems in protein structure prediction to contribute to folding. OPUS-Fold is based on torsion-angle sampling. After each sampling step, it reconstructs the structure and estimates the model quality with an energy function that is formed by combining many different constraining terms designed either by ourselves or by others in literature. OPUS-Fold balances the accuracy and the efficiency, delivers good results in a short time, and leaves more space for including the results of other subproblem methods. Moreover, OPUS-Fold also contains a fast side-chain modeling method OPUS-Rota2 (J. Chem. Theory Comput. 2019, 15 (9), 5154-5160), which enables a speedy construction of all-atom atomic models during the folding process that allows the usage of all-atom-required subproblem methods. In summary, OPUS-Fold provides a protein folding platform for incorporating the results from various subproblem methods, including those containing non-differentiable information such as partial experimental data.

## TBD
