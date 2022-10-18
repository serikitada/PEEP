# Metrics 
## `--metrics` (metrics to be used in filtering, ranking or to be shown in the output)

### On-target activity
* `DeepPE` 
* `DeepSpCas9`
* `CRISPRscan`
* `RuleSet1`
### Off-target activity
* `CFDscore`
* `MITscore`
* `mismatch_hit`

### Miscellaneous
* `chromatin_contact`
accpets bedpe. or built-in library
* `PolyT`
* `Spacer_GC`
* `PBS_GC`
* `PBS_RTT_GCrich`


## `--filter` 
### On-target activity
* `DeepPE` 
* `DeepSpCas9`
* `CRISPRscan`
* `RuleSet1`
### Off-target activity
* `CFDscore`
* `MITscore`
* `mismatch_hit`

## `--rank_each`
### On-target activity
* `DeepPE` 
* `DeepSpCas9`
* `CRISPRscan`
* `RuleSet1`
### Off-target activity
* `CFDscore`
* `MITscore`

## `--rank_pair`
"sum"
"product"
"min"

##flags
Poly T ( ≧ 4 consecutive Ts) flag (RNA polymerase III terminator)
Spacer GC flag (≦25%, ≧75%)
PBS GC flag (≦30%, ≧60%)
3’ extension (PBS+RTT) GC rich ( ≧ 4 GC counts out of 5bp) flag*

