# Fixed unrooted quartet: JC branch-length probe

`(T1:0.1,T2:0.2,(T3:0.3,T4:0.4):0.5);`

200 independent MSAs × 5000 sites per backend. SE is calculated across MSA replicate means.

| Backend | Tip pair | Tree distance | JC expectation | Observed mismatch | SE |
| --- | --- | ---: | ---: | ---: | ---: |
| bridge | T1_T2 | 0.3 | 0.24726 | 0.24761 | 0.00042 |
| bridge | T1_T3 | 0.9 | 0.52410 | 0.52341 | 0.00049 |
| bridge | T1_T4 | 1.0 | 0.55230 | 0.55220 | 0.00049 |
| bridge | T2_T3 | 1.0 | 0.55230 | 0.55229 | 0.00049 |
| bridge | T2_T4 | 1.1 | 0.57698 | 0.57711 | 0.00052 |
| bridge | T3_T4 | 0.7 | 0.45507 | 0.45410 | 0.00050 |
| standalone_240 | T1_T2 | 0.3 | 0.24726 | 0.24702 | 0.00041 |
| standalone_240 | T1_T3 | 0.9 | 0.52410 | 0.52443 | 0.00051 |
| standalone_240 | T1_T4 | 1.0 | 0.55230 | 0.55307 | 0.00050 |
| standalone_240 | T2_T3 | 1.0 | 0.55230 | 0.55234 | 0.00050 |
| standalone_240 | T2_T4 | 1.1 | 0.57698 | 0.57730 | 0.00048 |
| standalone_240 | T3_T4 | 0.7 | 0.45507 | 0.45565 | 0.00048 |
| standalone_2226 | T1_T2 | 0.3 | 0.24726 | 0.17574 | 0.00040 |
| standalone_2226 | T1_T3 | 0.9 | 0.52410 | 0.49139 | 0.00051 |
| standalone_2226 | T1_T4 | 1.0 | 0.55230 | 0.52410 | 0.00047 |
| standalone_2226 | T2_T3 | 1.0 | 0.55230 | 0.55144 | 0.00050 |
| standalone_2226 | T2_T4 | 1.1 | 0.57698 | 0.57701 | 0.00046 |
| standalone_2226 | T3_T4 | 0.7 | 0.45507 | 0.45494 | 0.00051 |

The bridge and standalone 2.4.0 track the JC expectations. In the installed
2.2.2.6 build, pairs involving T1 instead track the predictions after subtracting
T1’s 0.1 terminal edge: 0.17555, 0.49188 and 0.52410. Other pairs track the original
distances. This is consistent with that build ignoring the first terminal edge
when rooting this unrooted input; the experiment does not establish whether the
cause is the version, build configuration, or platform.

Reproduce with `python scripts/probe_alisim_branches.py --output NEW_DIRECTORY` in the validation environment. This is a descriptive diagnostic; it is not included in the main report’s family of KS checks.
