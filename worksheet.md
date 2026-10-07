# Worksheet: sleep deprivation and gene expression in mouse prefrontal cortex

*Omics for Neuroscience, computer practical*

Name / ID: ______________________________

Work through the notebook from top to bottom. In each section, make the prediction before
changing any settings, then answer the questions. Write short answers here or in the
notebook's answer boxes, and download them at the end.

---

## 1. The experiment and the data

Prediction: after 5 h of sleep loss, *Arc* expression will be ☐ higher ☐ lower ☐ unchanged

1. Was your prediction correct? By roughly how many times did *Arc* change?

2. Choose *Shank3* with raw counts, then with CPM. What changed, and why? (Hint: the
   library-size chart.)

3. How much do the 5 mice within one group differ for *Fos*? Why do we need several mice?

---

## 2. Testing thousands of genes

Prediction: 10,000 genes, none truly change, raw *p* < 0.05. How many will be "significant"?
☐ ~0 ☐ ~50 ☐ ~500 ☐ ~5,000

4. With 0% true changes, what happened when you switched from raw p-values to BH-adjusted
   p-values?

5. With 10% true changes and BH at 0.01 / 0.05 / 0.2, how close was the false discovery
   proportion to the threshold?

6. Real data, Shank3 vs WT: how many genes at raw *p* < 0.05, how many expected by chance, how
   many after BH? What did the shuffle show?

7. Describe the p-value histogram for sleep deprivation. What does the spike near 0 mean?

---

## 3. Differential expression

8. Raising the minimum |log2 fold change| to 1: how many genes remain, and which marker genes
   drop out? Is a small but very reliable change unimportant?

9. Find a gene with a big fold change that is *not* significant. Why isn't it?

10. *Account for day*: how many extra sleep-loss genes? Which genotype genes appear?

11. Is *Shank3* itself differentially expressed in the mutants? Why might it not be?

12. How does the number of discoveries change with the number of mice per group?

---

## 4. Gene Ontology enrichment

Prediction: which biological processes will sleep loss change?
______________________________

13. Top GO term for up-regulated genes with |log2FC| ≥ 0, and with |log2FC| ≥ 1. Why do they
    differ? Which result would you report?

14. "Response to glucocorticoid": is it caused by sleep loss or by the stress of being kept
    awake? How could an experiment tell?

15. What happened with the "all genes" background? Why is it wrong here?

16. How many terms did random gene lists give? What does that tell you?

---

## 5. Principal component analysis

Prediction: PC1 will separate ☐ genotype ☐ condition ☐ experiment day ☐ none of these

17. What does PC1 separate? What drives PC2/PC3?

18. Do WT and Shank3 mice separate? Does that mean the mutation has no effect?

19. What changed without the log transformation, and why?

20. In Section 5b, why do the random groups separate when genes are selected using them?

---

## 6. Report (one paragraph)

What does 5 h of sleep loss do to the prefrontal cortex transcriptome? How confident are you,
and what are the main caveats? Does *Shank3* matter?
