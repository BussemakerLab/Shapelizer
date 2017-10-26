# Shapelizer
---

Current purposes of this repository is:
* To distribute the methods for analyzing DNA shape readout developed in Rube et al. 2017:
* To illustrate how these methods were applied to generate the data shown in figures 2,3,5 of the paper.

---

The files in the repository are:
* `package/computeMeanProfile.py`        - Computes the mean shape profile of given a list of sequences
* `package/package/sampleSequences.py`   - Samples random sequences, computes the affinity for each sequence, and bins into  affinity bins
* `package/kMerLinearRegression.py`      - Performs linear regression on a k-mer table and outputs r, cross-validated r, regression coefficients, or the predicted values
* `package/linearModelToKMer.py`         - Creates a k-mer table by evaluating a linear model
* `package/shapeProjection.py`           - Performs shape projection given a mono+di binding model and a mono+di shape model
* `package/shapeLibrary.py`              - Various functions

* `scripts/fig2.sh`                      - Generates data for Figure 2 (run make fig2)
* `scripts/fig3.sh`                      - Generates data for Figure 3 (run make fig3)
* `scripts/fig5.sh`                      - Generates data for Figure 5 (run make fig5)

* `input/bindingModels`                  - TF binding models used in Figure 2 and 5
* `input/highAffinitySequences`          - High-affinity sequences used for Supplemental figure 1
* `input/shapeTables`                    - Shape tables
