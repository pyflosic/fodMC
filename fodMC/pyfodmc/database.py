def write_database(path='./'):
    data = '''############################################
# All average radii for the atoms (in bohr). 
# Structured as follows:
# 
# Element_symbol	Element_number		Pseudo_potential_electrons	No_of_core_electrons		Covalent_radius		Comment
# number_of_UP_shells   number_of_DN_shells	always at least 1. The number of points on this shell can be zero (if one spin channel shall be neglected)
# radius_UP_1   	number_of_points_UP_1
# radius_UP_2   	number_of_points_UP_2
# ....
# radius_DN_1   	number_of_points_DN_1
# ...
#############################################
##################
## ALL-ELECTRON ##
##################
H	1	0	0       0.605	Hydrogen atom. 	1 UP, 0 DN     	
1 1
0.300   1
0.300   0
He      2	0	0	1.757	Helium atom.   	1 UP, 1 DN	
1 1
0.300   1
0.300   1
Li      3	0	2	2.324	Lithium atom.  	2 UP, 1 DN	ADD aritifical outer (empty) shell to assign core and valence correctly
2 2
0.138   1
2.887   1
0.101   1
2.887   0
Be      4	0	2	1.701	Beryllium atom. 2 UP, 2 DN	
2 2
0.078   1
2.038   1
0.078   1
2.038   1
B       5	0	2	1.549	Boron atom.	3 UP, 2 DN
2 2
0.035   1
2.013   2
0.023   1
2.010   1
C       6	0	2	1.455	Carbon atom.	4 UP, 2 DN
2 2
0.083   1
2.001   3
0.069   1
2.003   1
C_3p    6       0       2	1.455	Carbon atom.    2 UP, 1 DN
2 2
0.083   1
2.001   1
0.069   1
2.003   0
C_HS    6       0       2	1.455	Carbon atom.    5 UP, 1 DN
2 2
0.083   1
2.001   4
0.069   1
2.003   0
N       7	0	2	1.417	Nitrogen atom.	5 UP, 2 DN
2 2
0.001	1
1.559	4
0.053	1
1.499	1
O	8	0	2	1.379	Oxygen atom.	5 UP, 3 DN
2 2
0.003	1
1.358	4
0.046	1
1.511	2
O_UP    8       0	2	1.379	Oxygen atom.    4 UP, 4 DN
2 2
0.010   1
1.450   3
0.010   1
1.450   3
F	9	0	2	1.361	Fluorine atom.	5 UP, 4 DN
2 2
0.002	1
1.204	4
0.004	1
1.406	3
Ne	10	0	2	1.342	Neon atom.	5 UP, 5 DN
2 2
0.001	1
1.079	4
0.001	1
1.079	4
Na	11	0	10	2.910	Sodium atom.	6 UP, 5 DN	ADD aritifical outer (empty) shell to assign core and valence correctly
3 3
0.002	1
0.882	4
4.829	1
0.002	1
0.916	4
4.829   0
Mg	12	0	10	2.570	Magnesium atom.	6 UP, 6 DN
3 3
0.001	1
0.711	4
3.920	1
0.001	1
0.711	4
3.920	1
Al	13	0	10	2.229	Aluminum atom.	7 UP, 6 DN
3 3
0.001	1
0.633	4
2.279	2
0.001	1
0.636	4
2.950	1	
Si	14	0	10	2.098	Silicon atom.	8 UP, 6 DN
3 3
0.001	1
0.558	4
2.051	3
0.001	1
0.559	4
2.555	1
P	15	0	10	2.003	Phosphorus atom. 9 UP, 6 DN
3 3
0.001	1
0.496	4
1.819	4
0.001	1
0.502	4
2.425	1
S	16	0	10	1.928	Sulfur atom.	9 UP, 7 DN
3 3
0.001	1
0.463	4
1.612	4
0.001	1
0.468	4
1.609	2
Cl	17	0	10	1.871	Chlorine atom.	9 UP, 8 DN
3 3
0.001	1
0.416	4
1.455	4
0.001	1
0.419	4
1.480	3
Ar	18	0	10	1.852	Argon atom.	9 UP, 9 DN
3 3
0.001	1
0.387	4
1.330	4
0.001	1
0.387	4
1.330	4
K	19	0	18	3.836	Potassium atom. 10 UP, 9 DN	ADD aritifical outer (empty) shell to assign core and valence correctly
4 4
0.001 	1
0.355	4
1.225	4
6.193	1
0.001	1
0.355	4
1.203	4
6.193	0
Ca	20	0	18	3.288	Calcium atom.	10 UP, 10 DN
4 4
0.001	1
0.323	4
1.285	4
3.957	1
0.001	1
0.323	4
1.284	4
3.999	1
Sc_4s1  21	0	18	2.721	Scandium atom.  12 UP, 9 DN    Standard for Sc is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d2 4s1 configuration
4 4
0.001   1
0.299   4
1.129   6
4.754   1
0.001   1
0.305   4
1.124   4
4.024   0
Sc_4s1_DN  21   0       18	2.721	Scandium atom.  9 UP, 12 DN    Standard for Sc is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d2 4s1 configuration
4 4
0.001   1
0.299   4
1.129   4
4.024   0
0.001   1
0.305   4
1.124   6
4.754   1
Sc_4s2  21      0       18	2.721	Scandium atom.  11 UP, 10 DN    Standard for Sc is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d1 4s2 configuration
4 4
0.001   1
0.299   4
1.129   5
4.754   1
0.001   1
0.305   4
1.124   4
4.024   1
Sc_4s2_DN  21   0       18	2.721	Scandium atom.  10 UP, 11 DN    Standard for Sc is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d1 4s2 configuration
4 4
0.001   1
0.299   4
1.124   4
4.024   1
0.001   1
0.305   4
1.129   5
4.754   1
Sc_3d4s	21	0	18	2.721	Scandium atom.	11 UP, 10 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.  Outer shell at arbitrary distance for now.
4 4
0.001   1
0.299   4
1.129   4
3.500   2
0.001   1
0.305   4
1.124   4
3.500   1
Ti_4s1  22	0	18	2.494	Titanium atom.  13 UP, 9 DN    	Standard for Ti is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d3 4s1 configuration
4 4
0.001   1
0.285   4
1.010   7
4.340   1
0.001   1
0.289   4
0.999   4
4.393   0
Ti_4s1_DN  22   0       18	2.494	Titanium atom.  9 UP, 13 DN     Standard for Ti is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d3 4s1 configuration
4 4
0.001   1
0.285   4
0.999   4
4.393   0
0.001   1
0.289   4
1.010   7
4.340   1
Ti_4s2  22      0       18	2.494	Titanium atom.  12 UP, 10 DN    Standard for Ti is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d2 4s2 configuration
4 4
0.001   1
0.285   4
1.010   6
4.340   1
0.001   1
0.289   4
0.999   4
4.393   1
Ti_4s2_DN  22   0       18	2.494	Titanium atom.  10 UP, 12 DN    Standard for Ti is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d2 4s2 configuration
4 4
0.001   1
0.285   4
0.999   4
4.393   1
0.001   1
0.289   4
1.010   6
4.340   1
Ti_3d4s	22	0	18	2.494	Titanium atom.	12 UP, 10 DN.	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.285   4
1.010   4
3.500   3
0.001   1
0.289   4
0.999   4
3.500   1
V_4s1   23	0	18	2.305	Vanadium atom.  14 UP, 9 DN    	Standard for V is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d4 4s1 configuration
4 4
0.001   1
0.261   4
1.025   8
4.538   1
0.001   1
0.271   4
1.026   4
4.019   0
V_4s1_DN   23   0       18	2.305	Vanadium atom.  9 UP, 14 DN     Standard for V is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d4 4s1 configuration
4 4
0.001   1
0.261   4
1.026   4
4.019   0
0.001   1
0.271   4
1.025   8
4.538   1
V_4s2   23      0       18	2.305	Vanadium atom.  13 UP, 10 DN    Standard for V is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d3 4s2 configuration
4 4
0.001   1
0.261   4
1.025   7
4.538   1
0.001   1
0.271   4
1.026   4
4.019   1
V_4s2_DN   23   0       18	2.305	Vanadium atom.  10 UP, 13 DN    Standard for V is 1s, 2s2p, 3s3p3d and 4s as different shells.  3d3 4s2 configuration
4 4
0.001   1
0.261   4
1.026   4
4.019   1
0.001   1
0.271   4
1.025   7
4.538   1
V_3d4s  23	0	18	2.305	Vanadium atom.  13 UP, 10 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells. 
4 4
0.001   1
0.261   4
1.025   4
4.538   4
0.001   1
0.271   4
1.026   4
4.019   1
Cr_4s1	24	0	18	2.229	Chromium atom.	15 UP, 9 DN.	Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.   3d5 4s1 configuration
4 4
0.001   1
0.253   4
1.013   9
3.978   1
0.001   1
0.268   4
1.044   4
4.000   0
Cr_4s1_DN  24   0       18	2.229	Chromium atom.  9 UP, 15 DN.    Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.   3d5 4s1 configuration
4 4
0.001   1
0.253   4
1.044   4
4.000   0
0.001   1
0.268   4
1.013   9
3.978   1
Cr_4s2  24      0       18	2.229	Chromium atom.  14 UP, 10 DN.    Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.   3d4 4s2 configuration
4 4
0.001   1
0.253   4
1.013   8
3.978   1
0.001   1
0.268   4
1.044   4
4.000   1
Cr_4s2_DN  24   0       18	2.229	Chromium atom.  10 UP, 14 DN.    Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.   3d4 4s2 configuration
4 4
0.001   1
0.253   4
1.044   4
4.000   1
0.001   1
0.268   4
1.013   8
3.978   1
Cr_3d4s 24      0       18	2.229	Chromium atom.  15 UP, 9 DN.    other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.253   4
1.013   4
3.978   6
0.001   1
0.268   4
1.044   4
4.000   0
Cr_6bond 24      0       18	2.229	Chromium atom.  15 UP, 9 DN.    other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.250   4
1.013   4
3.978   6
0.001   1
0.250   4
1.044   4
4.000   0
Cr_test1  24     0	18	2.229	Chromium atom.  14 UP, 10 DN.  Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.
4 4
0.001   1
0.253   4
0.450   6
1.013   4
0.001   1
0.268   4
1.044   4
4.000   0
Cr_test2  24	0       18	2.229	Chromium atom.  10 UP, 14 DN. 	Standard for Cr is 1s, 2s2p, 3s3p3d and 4s as different shells.   3d5 4s1 configuration
4 4
0.001   1
0.253   4
1.044   4
4.000   0
0.001   1
0.268   4
0.450   6
1.031   4
Mn_4s1  25	0	18	2.211	Manganese atom. 15 UP, 10 DN    Standard for Mn is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s1 configuration
4 4
0.001   1
0.236   4
0.904   9
6.067   1
0.001   1
0.249   4
0.905   5
3.945   0
Mn_4s1_DN  25   0       18	2.211	Manganese atom. 10 UP, 15 DN    Standard for Mn is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s1 configuration
4 4
0.001   1
0.236   4
0.905   5
3.945   0
0.001   1
0.249   4
0.904   9
6.067   1
Mn_4s2  25      0       18	2.211	Manganese atom. 15 UP, 10 DN    Standard for Mn is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d5 4s2 configuration
4 4
0.001   1
0.236   4
0.904   9
6.067   1
0.001   1
0.249   4
0.905   4
3.945   1
Mn_4s2_DN  25   0       18	2.211	Manganese atom. 10 UP, 15 DN    Standard for Mn is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d5 4s2 configuration
4 4
0.001   1
0.236   4
0.905   4
3.945   1
0.001   1
0.249   4
0.904   9
6.067   1
Mn_3d4s 25	0	18	2.211	Manganese atom. 15 UP, 10 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.236   4
0.904   4
6.067   6
0.001   1
0.249   4
0.905   4
3.945   1
Fe_4s1  26	0	18	2.211	Iron atom.      15 UP, 11 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d7 4s1 configuration
4 4
0.001   1
0.226   4
0.849   9
4.660   1
0.001   1
0.240   4
0.882   6
3.238   0
Fe_4s1_DN  26   0       18	2.211	Iron atom.      11 UP, 15 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d7 4s1 configuration
4 4
0.001   1
0.226   4
0.882   6
3.238   0
0.001   1
0.240   4
0.849   9
4.660   1
Fe_4s2  26      0       18	2.211	Iron atom.      15 UP, 11 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   9
4.660   1
0.001   1
0.240   4
0.882   5
3.238   1
Fe_4s2_DN  26      0    18	2.211	Iron atom.      11 UP, 15 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.882   5
3.238   1
0.001   1
0.240   4
0.849   9
4.660   1
Fe_4s2_HS3  26     0   18	2.211	Iron atom.      14 UP, 12 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   8
4.660   1
0.001   1
0.240   4
0.882   6
3.238   1
Fe_4s2_LS1  26     0   18	2.211	Iron atom.      13 UP, 13 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   7
4.660   1
0.001   1
0.240   4
0.882   7
3.238   1
Fe_4s2_LS1_2+  26     0   18	2.211		Iron ion.      13 UP, 13 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
3 3
0.001   1
0.226   4
1.249   7
0.001   1
0.240   4
1.282   7
Fe_4s2_HS5_2+  26     0   18	2.211		Iron ion.      15 UP, 11 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
3 3
0.001   1
0.226   4
1.249   9
0.001   1
0.240   4
1.282   5
Fe_4s2_6bonds_HS5  26  0  18	2.211     	Iron atom.   15 UP, 11 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   7
4.660   3
0.001   1
0.240   4
0.882   3
3.238   3
Fe_4s2_6bonds_HS3  26  0   18	2.211    	Iron atom.   14 UP, 12 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.   
4 4
0.001   1
0.226   4
0.849   6
4.660   3
0.001   1
0.240   4
0.882   4
3.238   3
Fe_4s2_6bonds_LS1  26  0   18	2.211    	Iron atom.   13 UP, 13 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    LS core
4 4
0.001   1
0.226   4
0.849   5
4.660   3
0.001   1
0.240   4
0.882   5
3.238   3
Fe_3d4s 26	0	18	2.211		Iron atom.      15 UP, 11 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.226   4
0.849   4
4.660   6
0.001   1
0.240   4
0.882   4
3.238   2
Fe_4s2_III  26   0	18	2.211       	Iron ion.      15 UP, 10 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   9
4.660   1
0.001   1
0.240   4
0.882   4
3.238   1
Fe_4s2_III_DN  26   0	18	2.211       	Iron ion.      10 UP, 15 DN    Standard for Fe is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d6 4s2 configuration
4 4
0.001   1
0.226   4
0.849   4
4.660   1
0.001   1
0.240   4
0.882   9
3.238   1
Fe_3s3p3d4s  	26  0  18	2.211     	Iron atom.      15 UP, 11 DN    1s, 2s2p, 3s3p3d4s as different shells.    3d6 4s2 configuration
3 3
0.001   1
0.226   4
0.849   10
0.001   1
0.240   4
0.882   6
Fe_test1 26      0	18	2.211       	Iron atom.      15 UP, 11 DN    other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.   for sebbi
4 4
0.001   1
0.226   4
0.849   4
1.300   6
0.001   1
0.240   4
0.882   4
1.300   2
Co_4s1  27	0	18	2.192	Cobalt atom.    15 UP, 12 DN    Standard for Co is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d8 4s1 configuration
4 4
0.001   1
0.216   4
0.829   9
4.360   1
0.001   1
0.235   4
0.822   7
3.517   0
Co_4s1_DN  27   0       18	2.192	Cobalt atom.    12 UP, 15 DN    Standard for Co is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d8 4s1 configuration
4 4
0.001   1
0.216   4
0.822   7
3.517   0
0.001   1
0.235   4
0.829   9
4.360   1
Co_4s2  27      0       18	2.192	Cobalt atom.    15 UP, 12 DN    Standard for Co is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d7 4s2 configuration
4 4
0.001   1
0.216   4
0.829   9
4.360   1
0.001   1
0.235   4
0.822   6
3.517   1
Co_4s2_DN  27   0       18	2.192	Cobalt atom.    12 UP, 15 DN    Standard for Co is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d7 4s2 configuration
4 4
0.001   1
0.216   4
0.822   6
3.517   1
0.001   1
0.235   4
0.829   9
4.360   1
Co_3d4s 27	0	18	2.192	Cobalt atom.    15 UP, 12 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.216   4
0.829   4
4.360   6
0.001   1
0.235   4
0.822   4
3.517   3
Ni_4s1  28	0	18	2.173	Nickel atom.	15 UP, 13 DN    Standard for Ni is 1s, 2s2p, 3s3p3d and 4s as different shells.     3d9 4s1 configuration
4 4
0.001   1
0.201   4
0.825   9
3.722   1
0.001   1
0.216   4
0.815   8
4.924   0
Ni_4s1_DN  28   0       18	2.173	Nickel atom.    13 UP, 15 DN    Standard for Ni is 1s, 2s2p, 3s3p3d and 4s as different shells.     3d9 4s1 configuration
4 4
0.001   1
0.201   4
0.815   8
4.924   0
0.001   1
0.216   4
0.825   9
3.722   1
4.924   0
Ni_4s2  28      0       18	2.173	Nickel atom.    15 UP, 13 DN    Standard for Ni is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d8 4s2 configuration
4 4
0.001   1
0.201   4
0.825   9
3.722   1
0.001   1
0.216   4
0.815   7
4.924   1
Ni_4s2_DN  28   0       18	2.173	Nickel atom.    13 UP, 15 DN    Standard for Ni is 1s, 2s2p, 3s3p3d and 4s as different shells.    3d8 4s2 configuration
4 4
0.001   1
0.201   4
0.815   7
4.924   1
0.001   1
0.216   4
0.825   9
3.722   1
Ni_3d4s 28	0	18	2.173	Nickel atom.	15 UP, 13 DN	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.201   4
0.825   4
3.722   6
0.001   1
0.216   4
0.815   4
4.924   4
Cu_4s1	29	0	18	2.211	Copper atom.	15 UP, 14 DN	Standard for Cu is 1s, 2s2p, 3s3p3d and 4s as different shells.      3d10 4s1  configuration
4 4
0.001	1
0.197	4
1.133	9
6.041	1
0.001	1
0.190	4
1.128	9
4.000	0
Cu_4s1_DN  29   0       18	2.211	Copper atom.    14 UP, 15 DN    Standard for Cu is 1s, 2s2p, 3s3p3d and 4s as different shells.      3d10 4s1  configuration
4 4
0.001   1
0.197   4
1.128   9
4.000   0
0.001   1
0.190   4
1.133   9
6.041   1
Cu_4s2  29	0	18	2.211	Copper atom.	15 UP, 14 DN    Standard for Cu is 1s, 2s2p, 3s3p3d and 4s as different shells.      3d9 4s2  configuration
4 4
0.001   1
0.197   4
1.133   9
6.041   1
0.001   1
0.190   4
1.128   8
5.000   1
Cu_4s2_DN  29   0       18	2.211	Copper atom.    14 UP, 15 DN    Standard for Cu is 1s, 2s2p, 3s3p3d and 4s as different shells.      3d9 4s2  configuration
4 4
0.001   1
0.197   4
1.128   8
5.000   1
0.001   1
0.190   4
1.133   9
6.041   1
Cu_3d4s 29	0	18	2.211	Copper atom.	15 UP, 14 DN    other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.197   4
1.133   4
6.041   6
0.001   1
0.190   4
1.128   4
5.000   5
Cu_test  29      0       18	2.211	Copper atom.    15 UP, 14 DN    Standard for Cu is 1s, 2s2p, 3s3p3d and 4s as different shells.      3d10 4s1  configuration
4 4
0.001   1
0.197   4
0.900   9
6.041   1
0.001   1
0.190   4
0.900   9
4.000   0
Zn	30	0	18	2.362	Zinc atom	15 UP, 15 DN	Standard for Zn is 1s, 2s2p, 3s3p3d and 4s as different shells.
4 4
0.001	1
0.182	4
0.755	9
4.745	1
0.001   1
0.182   4
0.755   9
4.745   1
Zn_3d4s	30	0	18	2.362	Zinc atom.	15 UP, 15 DN 	other guess. Taking 1s, 2s2p, 3s3p and 3d4s as different shells.
4 4
0.001   1
0.182   4
0.755   4
4.745   6
0.001   1
0.182   4
0.755   4
4.745   6
Ga	31	0	28	2.381	Gallium atom.	16 UP, 15 DN.	
4 4
0.001	1
0.180	4
0.706	9
4.165	2
0.001	1
0.180	4
0.712	9
4.325	1
Ge	32	0	28	2.305	Germanium atom.	17 UP, 15 DN	
4 4
0.001	1
0.174	4
0.670	9
3.334	3
0.001	1
0.173	4
0.676	9
3.296	1
As	33	0	28	2.268	Arsen atom.	18 UP, 15 DN	
4 4
0.001	1
0.165	4
0.623	9
3.202	4
0.001	1
0.166	4
0.628	9
2.953	1
Se	34	0	28	2.192	Selenium atom.	18 UP, 16 DN	
4 4
0.001 	1
0.159	4
0.582	9
3.192	4
0.001	1
0.159	4
0.590	9
2.910	2
Br	35	0	28	2.154	Bromium atom.	18 UP, 17 DN	
4 4
0.001	1
0.153	4
0.556	9
2.647	4
0.001	1
0.154	4
0.564	9
2.634	3
Kr	36	0	28	2.116	Krypton atom.	18 UP, 18 DN	
4 4
0.001	1
0.147	4
0.523	9
2.475	4
0.001   1
0.147   4
0.523   9
2.475   4
Kr_test 36	0	28	2.116	Krypton atom.   18 UP, 18 DN
5 5
0.000   1
0.147   4
0.523   4
1.500   5
2.475   4
0.000   1
0.147   4
0.523   4
1.500   5
2.475   4
Zr      40      0       36      0.10  -  (2.797)_org   Zirconium atom.   21 UP, 19 DN  -- don't give it any bonds
5 5
0.001   1
0.147   4
0.523   9
1.475   6
3.500   1
0.001   1
0.147   4
0.523   9
1.475   4
3.500   1
Zr_LS   40      0       36      0.10  -  (2.797)_org   Zirconium atom.   20 UP, 20 DN  -- don't give it any bonds
5 5
0.001   1
0.147   4
0.523   9
1.475   5
3.500   1
0.001   1
0.147   4
0.523   9
1.475   5
3.500   1
Xe      54     0	46	2.476	Xenon atom.	27 UP, 27 DN
5 5
0.000	1
0.095	4
0.273	9
0.848	9
2.678	4
0.000	1
0.095	4
0.273	9
0.848	9
2.678	4
Rn      86     0        78	10.0	Radon atom.     43 UP, 43 DN
6 6
0.000   1
0.095   4
0.273   9
0.848   9
1.100   16
2.678   4
0.000   1
0.095   4
0.273   9
0.848   9
1.100   16
2.678   4
###########################
## LC ECPs - large core ##
###########################
H_ECP_LC       	1	0	0	0.605 	Hydrogen atom.  1s1   		ECP guess	1 UP 0 DN
1 1
0.300   1
0.300   0
He_ECP_LC      	2 	0	0	1.757	Helium atom.    1s2   		ECP guess	1 UP 1 DN
1 1
0.300   1
0.300   1
Li_ECP_LC      	3	2	2	2.324	Lithium atom.   2s1		ECP guess      	1 UP 0 DN	ADD aritifical outer (empty) shell to assign core and valence correctly
1 1
2.887   1
2.887   0
Be_ECP_LC      	4 	2	2	1.701	Beryllium atom. 2s2		ECP guess	1 UP 1 DN
1 1
2.038   1
2.038   1
B_ECP_LC	5	2 	2	1.549	Boron atom.     2s2_2p1		ECP guess	2 UP 1 DN
1 1
2.013   2
2.010   1
C_ECP_LC       	6	2 	2	1.455	Carbon atom.    2s2_2p2		ECP guess	3 UP 1 DN
1 1
2.001   3
2.003   1
N_ECP_LC       	7 	2	2	1.417	Nitrogen atom.  2s2_2p3		ECP guess	4 UP 1 DN
1 1
1.559   4
1.499   1
O_ECP_LC       	8 	2	2	1.379	Oxygen atom.    2s2_2p4		ECP guess	4 UP 2 DN
1 1
1.358   4
1.511   2
F_ECP_LC       	9	2	2	1.361	Fluorine atom.  2s2_2p5		ECP guess	4 UP 3 DN
1 1
1.204   4
1.406   3
Ne_ECP_LC      	10	2 	2	1.342	Neon atom.     2s2_2p6		ECP guess	4 UP 4 DN
1 1
1.079   4
1.079   4
Na_ECP_LC      	11	10	10	2.910	Sodium atom.   3s1		ECP guess      	1 UP 0 DN	ADD aritifical outer (empty) shell to assign core and valence correctly
1 1
4.829   1
4.829   0
Mg_ECP_LC      	12	10 	10	2.570	Magnesium atom. 3s2		ECP guess	1 UP 1 DN
1 1
3.920   1
3.920   1
Al_ECP_LC	13	10 	10	2.229	Aluminum atom. 3s2_3p1		ECP guess	2 UP 1 DN
1 1
2.279   2
2.950   1
Si_ECP_LC      	14	10 	10	2.098	Silicon atom.  3s2_3p2		ECP guess	3 UP 1 DN
1 1
2.051   3
2.555   1
P_ECP_LC       	15	10 	10	2.003	Phosphorus atom. 3s2_3p3	ECP guess	4 UP 1 DN
1 1
1.819   4
2.425   1
S_ECP_LC       	16	10 	10	1.928	Sulfur atom.   3s2_3p4		ECP guess	4 UP 2 DN
1 1
1.612   4
1.609   2
Cl_ECP_LC      	17	10 	10	1.871	Chlorine atom. 3s2_3p5		ECP guess	4 UP 3 DN
1 1
1.455   4
1.480   3
Ar_ECP_LC	18	10 	10	1.852	Argon atom.    3s2_3p6		ECP guess	4 UP 4 DN
1 1
1.330   4
1.330   4
K_ECP_LC       	19	18 	18	3.386	Potassium atom. 4s1		ECP guess	1 UP 0 DN 	ADD aritifical outer (empty) shell to assign core and valence correctly
1 1
6.193   1
6.193   0
Ca_ECP_LC	20	18 	18	3.288	Calcium atom.   4s2		ECP guess	1 UP 1 DN
1 1
3.957   1
3.999   1
Ga_ECP_LC      	31	28 	28	2.381	Gallium atom.   4s2_4p1		ECP guess	2 UP 1 DN
1 1
4.165   2
4.325   1
Ge_ECP_LC      	32 	28	28	2.305	Germanium atom  4s2_4p2		ECP guess	3 UP 1 DN
1 1
3.334   3
3.296   1
As_ECP_LC      	33	28 	28	2.268	Arsen atom.     4s2_4p3		ECP guess	4 UP 1 DN 
1 1	
3.202   4
2.953   1
Se_ECP_LC      	34	28 	28	2.192	Selenium atom.  4s2_4p4		ECP guess	4 UP 2 DN
1 1
3.192   4
2.910   2
Br_ECP_LC      	35	28 	28	2.154	Bromium atom.   4s2_4p5		ECP guess	4 UP 3 DN
1 1
2.647   4
2.634   3
Kr_ECP_LC      	36	28 	28	2.116	Krypton atom.   4s2_4p6		ECP guess	4 UP 4 DN
1 1
2.475   4
2.475   4
Rb_ECP_LC      	37	36 	36	Rubidium atom.  5s1    	ECP guess       1 UP 0 DN       ADD aritifical outer (empty) shell to assign core and valence. FROM here: values for radii are extrapolated from lighter elements
1 1
10.009  1
10.009  0
Sr_ECP_LC      	38	36 	36	Strontium atom. 5s2		ECP guess	1 UP 1 DN
1 1
6.424   1
6.424   1
In_ECP_LC	49	46 	46	Indium atom.    5s2_5p1		ECP guess	2 UP 1 DN
1 1
5.644   2
5.934   1
Sn_ECP_LC      	50 	46	46	Tin atom        5s2_5p2		ECP guess	3 UP 1 DN
1 1
4.274   3
4.195   1
Sb_ECP_LC      	51	46 	46	Antimony atom.  5s2_5p3		ECP guess	4 UP 1 DN
1 1
4.333   4
3.972   1
Te_ECP_LC      	52	46 	46	Tellurium atom. 5s2_5p4		ECP guess	4 UP 2 DN
1 1
4.454   4
3.869   2
I_ECP_LC       	53	46 	46	Iod atom.       5s2_5p5		ECP guess	4 UP 3 DN
1 1
3.641   4
3.476   3
Xe_ECP_LC      	54 	46	46	Xe atom.        5s2_5p6		ECP guess	4 UP 4 DN
1 1
3.593   4
3.593   4
Cs_ECP_LC      	55	54 	54	Caesium atom.   6s1		ECP guess	1 UP 0 DN	ADD aritifical outer (empty) shell to assign core and valence correctly.
1 1
13.729  1
13.729  0
Ba_ECP_LC      	56	54 	54	Barium atom.    6s2		ECP guess	1 UP 1 DN
1 1
8.583   1
8.583   1
Tl_ECP_LC      	81	78 	78	Thallium atom.  6s2_6p1		ECP guess	2 UP 1 DN
1 1
8.412   2
8.715   1
Pb_ECP_LC      	82 	78	78	Lead atom       6s2_6p2		ECP guess	3 UP 1 DN
1 1
6.035   3
5.739   1
Bi_ECP_LC      	83 	78	78	Bismuth atom.   6s2_6p3		ECP guess	4 UP 1 DN
1 1
6.429   4
5.617   1
Po_ECP_LC      	84	78 	78	Pollonium atom. 6s2_6p4		ECP guess	4 UP 2 DN
1 1
6.805   4
5.691   2
At_ECP_LC      	85 	78	78	Astate atom.    6s2_6p5		ECP guess	4 UP 3 DN
1 1
5.475   4
5.079   3
Rn_ECP_LC      	86 	78	78	Radon atom.     6s2_6p6		ECP guess	4 UP 4 DN
1 1
5.553   4
5.553   4
###################
# Small core ECPs #
###################
H_ECP_SC       	1	0 	0	0.605	Hydrogen atom.		Core = None					1 UP, 0 DN
1 1
0.300   1
0.300   0
He_ECP_SC      	2 	0	0	1.757	Helium atom.		Core = None					1 UP, 1 DN
1 1
0.300   1
0.300   1
Li_ECP_SC      	3	2 	2	2.324	Lithium atom.		Core = 1s		ECP guess small core   	1 UP, 0 DN      ADD aritifical outer (empty) shell to assign core and valence correctly
1 1
2.887   1
2.887   0
Be_ECP_SC      	4 	2	2	1.701	Beryllium atom. 	Core = 1s		ECP guess small core	1 UP, 1 DN
1 1
2.038   1
2.038   1
B_ECP_SC 	5	2 	2	1.549	Boron atom.		Core = 1s		ECP guess small core    2 UP, 1 DN
1 1
2.013   2
2.010   1
C_ECP_SC	6 	2	2	1.455	Carbon atom.		Core = 1s		ECP guess small core    3 UP, 1 DN
1 1
2.001   3
2.003   1
N_ECP_SC	7 	2	2	1.417	Nitrogen atom.		Core = 1s		ECP guess small core	4 UP, 1 DN
1 1
1.559   4
1.499   1
O_ECP_SC 	8	2 	2	1.379	Oxygen atom.		Core = 1s		ECP guess small core    4 UP, 2 DN
1 1
1.358   4
1.511   2
F_ECP_SC       	9	2 	2	1.361	Fluorine atom.		Core = 1s		ECP guess small core	4 UP, 3 DN
1 1
1.204   4
1.406   3
Ne_ECP_SC      	10	2 	2	1.342	Neon atom.		Core = 1s		ECP guess small core	4 UP, 4 DN
1 1
1.079   4
1.079   4
Na_ECP_SC	11 	10	10	2.910	Sodium atom.		Core = 1s, 2s2p		ECP guess small core    1 UP, 0 DN      ADD aritifical outer (empty) shell to assign core and valence correctly
1 1
4.829   1
4.829   0
Mg_ECP_SC      	12	10 	10	2.570	Magnesium atom. 	Core = 1s, 2s2p		ECP guess small core	1 UP, 1 DN
1 1
3.920   1
3.920   1
Al_ECP_SC      	13 	10	10	2.229	Aluminum atom.		Core = 1s, 2s2p		ECP guess small core	2 UP, 1 DN
1 1
2.279   2
2.950   1
Si_ECP_SC      	14	10 	10	2.098	Silicon atom.		Core = 1s, 2s2p		ECP guess small core	3 UP, 1 DN
1 1
2.051   3
2.555   1
P_ECP_SC	15	10 	10	2.003	Phosphorus atom. 	Core = 1s, 2s2p		ECP guess small core	4 UP, 1 DN
1 1
1.819   4
2.425   1
S_ECP_SC       	16 	10	10	1.928	Sulfur atom.		Core = 1s, 2s2p		ECP guess small core    4 UP, 2 DN
1 1
1.612   4
1.609   2
Cl_ECP_SC      	17	10 	10	1.871	Chlorine atom.		Core = 1s, 2s2p		ECP guess small core	4 UP, 3 DN
1 1
1.455   4
1.480   3
Ar_ECP_SC	18	10 	10	1.852	Argon atom.		Core = 1s, 2s2p		ECP guess small core	4 UP, 4 DN
1 1
1.330   4
1.330   4
K_ECP_SC       	19	10 	18	3.836	Potassium atom.		Core = 1s, 2s2p		ECP guess small core	5 UP, 4 DN     ADD aritifical outer (empty) shell to assign core and valence correctly
2 2
1.255   4
6.193   1
1.183   4
6.193   0
Ca_ECP_SC      	20	10 	18	3.288	Calcium atom.		Core = 1s, 2s2p		ECP guess small core	5 UP, 5 DN
2 2
1.285   4
3.957   1
1.284   4
3.999   1
Sc_ECP_SC_4s1	21	10 	18	2.721	Scandium atom.       	Core = 1s, 2s2p         ECP guess small core    7 UP, 4 DN      3d2 4s1 configuration
2 2
1.100   6
4.500   1
1.100   4
4.500   0
Sc_ECP_SC_4s1_DN   21   10      18	2.721	Scandium atom.          Core = 1s, 2s2p         ECP guess small core    4 UP, 7 DN      3d2 4s1 configuration
2 2
1.100   4
4.500   0
1.100   6
4.500   1
Sc_ECP_SC_4s2   21      10      18	2.721	Scandium atom.          Core = 1s, 2s2p         ECP guess small core    6 UP, 5 DN      3d1 4s2 configuration
2 2
1.100   5
4.500   1
1.100   4
4.500   1
Sc_ECP_SC_4s2_DN   21   10      18	2.721	Scandium atom.          Core = 1s, 2s2p         ECP guess small core    5 UP, 6 DN      3d1 4s2 configuration
2 2
1.100   4
4.500   1
1.100   5
4.500   1
Ti_ECP_SC_4s1	22	10	18	2.494	Titanium atom.		Core = 1s, 2s2p		ECP guess small core    8 UP, 4 DN      3d3 4s1 configuration
2 2
1.010   7
4.350   1
1.010   4
4.350   0
Ti_ECP_SC_4s1_DN   22   10      18	2.494	Titanium atom.          Core = 1s, 2s2p         ECP guess small core    4 UP, 8 DN      3d3 4s1 configuration
2 2
1.010   4
4.350   0
1.010   7
4.350   1
Ti_ECP_SC_4s2   22      10      18	2.494	Titanium atom.          Core = 1s, 2s2p         ECP guess small core    7 UP, 5 DN      3d2 4s2 configuration
2 2
1.010   6
4.350   1
1.010   4
4.350   1
Ti_ECP_SC_4s2_DN   22   10      18	2.494	Titanium atom.          Core = 1s, 2s2p         ECP guess small core    5 UP, 7 DN      3d2 4s2 configuration
2 2
1.010   4
4.350   1
1.010   6
4.350   1
V_ECP_SC_4s1	23	10 	18	2.305	Vanadium atom.		Core = 1s, 2s2p		ECP guess small core	9 UP, 4 DN	3d4 4s1 configuration
2 2
1.050   8
4.300   1
1.050   4
4.300   0
V_ECP_SC_4s1_DN    23   10      18	2.305	Vanadium atom.          Core = 1s, 2s2p         ECP guess small core    4 UP, 9 DN      3d4 4s1 configuration
2 2
1.050   4
4.300   0
1.050   8
4.300   1
V_ECP_SC_4s2    23      10      18	2.305	Vanadium atom.          Core = 1s, 2s2p         ECP guess small core    8 UP, 5 DN      3d3 4s2 configuration
2 2
1.050   7
4.300   1
1.050   4
4.300   1
V_ECP_SC_4s2_DN    23   10      18	2.305	Vanadium atom.          Core = 1s, 2s2p         ECP guess small core    5 UP, 8 DN      3d3 4s2 configuration
2 2
1.050   4
4.300   1
1.050   7
4.300   1
Cr_ECP_SC_4s1	24	10 	18	2.229	Chromium atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 4 DN.	3d5 4s1 configuration
2 2
1.100   9
4.000   1
1.100   4
4.000   0
Cr_ECP_SC_4s1_DN   24   10      18	2.229	Chromium atom.          Core = 1s, 2s2p         ECP guess small core    4 UP, 10 DN.    3d5 4s1 configuration
2 2
1.100   4
4.000   0
1.100   9
4.000   1
Cr_ECP_SC_4s2   24      10      18	2.229	Chromium atom.          Core = 1s, 2s2p         ECP guess small core    9 UP, 5 DN.     3d4 4s2 configuration
2 2
1.100   8
4.000   1
1.100   4
4.000   1
Cr_ECP_SC_4s2_DN   24   10      18	2.229	Chromium atom.          Core = 1s, 2s2p         ECP guess small core    5 UP, 9 DN.     3d4 4s2 configuration
2 2
1.100   4
4.000   1
1.100   8
4.000   1
Mn_ECP_SC_4s1	25	10 	18	2.211	Manganese atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 5 DN	3d6 4s1 configuration
2 2
0.900   9
4.500   1
0.900   5
4.500   0
Mn_ECP_SC_4s1_DN   25   10      18	2.211	Manganese atom.         Core = 1s, 2s2p         ECP guess small core    5 UP, 10 DN     3d6 4s1 configuration
2 2
0.900   5
4.500   0
0.900   9
4.500   1
Mn_ECP_SC_4s2   25      10      18	2.211	Manganese atom.         Core = 1s, 2s2p         ECP guess small core    10 UP, 5 DN     3d5 4s2 configuration
2 2
0.900   9
4.500   1
0.900   4
4.500   1
Mn_ECP_SC_4s2_DN   25   10      18	2.211	Manganese atom.         Core = 1s, 2s2p         ECP guess small core    5 UP, 10 DN     3d5 4s2 configuration
2 2
0.900   4
4.500   1
0.900   9
4.500   1
Mn_ECP_SC_d4    25    	10	18	2.211	rajendra. one d missing
2 2
0.900   8
4.500   1
0.900   4
4.500   1
Mn_ECP_SC_d3    25    	10	18	2.211	rajendra. two d missing
2 2
0.900   7
4.500   1
0.900   4
4.500   1
Fe_ECP_SC_4s1   26      10      18	2.211	Iron atom.              Core = 1s, 2s2p         ECP guess small core    10 UP, 6 DN     3d7 4s1 configuration
2 2
0.850   9
4.300   1
0.850   6
4.300   0
Fe_ECP_SC_4s1_DN   26   10      18	2.211	Iron atom.              Core = 1s, 2s2p         ECP guess small core    6 UP, 10 DN     3d7 4s1 configuration
2 2
0.850   6
4.300   0
0.850   9
4.300   1
Fe_ECP_SC_4s2	26 	10	18	2.211	Iron atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 6 DN	3d6 4s2 configuration
2 2 
0.850   9
4.300   1
0.850   5
4.300   1
Fe_ECP_SC_4s2_DN   26 	10	18	2.211	Iron atom.        	Core = 1s, 2s2p         ECP guess small core    6 UP, 10 DN     3d6 4s2 configuration
2 2
0.850   5
4.300   1
0.850   9
4.300   1
Fe_ECP_SC_4s2_HS3  26    10   	18	2.211	Iron atom.  		Core = 1s, 2s2p         ECP guess small core    9 UP, 7 DN
2 2
0.850   8
4.660   1
0.850   6
3.238   1
Fe_ECP_SC_4s2_LS1  26    10   	18	2.211	Iron atom.     		 Core = 1s, 2s2p         ECP guess small core    8 UP, 8 DN 
2 2
0.850   7
4.660   1
0.850   7
3.238   1
Co_ECP_SC_4s1   27      10      18	2.192	Cobalt atom.            Core = 1s, 2s2p         ECP guess small core    10 UP, 7 DN     3d8 4s1 configuration
2 2
0.810   9
4.000   1
0.810   7
4.000   0
Co_ECP_SC_4s1_DN   27   10      18	2.192	Cobalt atom.            Core = 1s, 2s2p         ECP guess small core    7 UP, 10 DN     3d8 4s1 configuration
2 2
0.810   7
4.000   0
0.810   9
4.000   1
Co_ECP_SC_4s2	27 	10	18	2.192	Cobalt atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 7 DN	3d7 4s2 configuration
2 2
0.810   9
4.000   1
0.810   6
4.000   1
Co_ECP_SC_4s2_DN   27   10      18	2.192	Cobalt atom.            Core = 1s, 2s2p         ECP guess small core    7 UP, 10 DN     3d7 4s2 configuration
2 2
0.810   6
4.000   1
0.810   9
4.000   1
Ni_ECP_SC_4s1   28      10      18	2.173	Nickel atom.            Core = 1s, 2s2p         ECP guess small core    10 UP, 8 DN     3d9 4s1 configuration
2 2
0.800   9
4.000   1
0.800   8
4.000   0
Ni_ECP_SC_4s1_DN   28   10      18	2.173	Nickel atom.            Core = 1s, 2s2p         ECP guess small core    8 UP, 10 DN     3d9 4s1 configuration
2 2
0.800   8
4.000   0
0.800   9
4.000   1
Ni_ECP_SC_4s2	28	10	 18	2.173	Nickel atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 8 DN	3d8 4s2 configuration
2 2
0.800   9
4.000   1
0.800   7
4.000   1
Ni_ECP_SC_4s2_DN   28   10       18	2.173	Nickel atom.           Core = 1s, 2s2p         ECP guess small core    8 UP, 10 DN     3d8 4s2 configuration
2 2
0.800   7
4.000   1
0.800   9
4.000   1
Cu_ECP_SC_4s1   29      10      18	2.211	Copper atom.            Core = 1s, 2s2p         ECP guess small core    10 UP, 9 DN     3d10 4s1 configuration
2 2
0.900   9
4.000   1
0.900   9
4.000   0
Cu_ECP_SC_4s1_DN   29   10      18	2.211	Copper atom.            Core = 1s, 2s2p         ECP guess small core    9 UP, 10 DN     3d10 4s1 configuration
2 2
0.900   9
4.000   0
0.900   9
4.000   1
Cu_ECP_SC_4s2	29	10	18	2.211	Copper atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 9 DN	3d9 4s2 configuration
2 2
0.900   9
4.000   1
0.900   8
4.000   1
Cu_ECP_SC_4s2_DN   29   10      18	2.211	Copper atom.            Core = 1s, 2s2p         ECP guess small core    9 UP, 10 DN     3d9 4s2 configuration
2 2
0.900   8
4.000   1
0.900   9
4.000   1
Zn_ECP_SC	30 	10	18	2.362	Zinc atom.		Core = 1s, 2s2p		ECP guess small core	10 UP, 10 DN	3d10 4s2 configuration
2 2
0.755   9
4.745   1
0.755   9
4.745   1
Ga_ECP_SC      	31	18 	28	2.381	Gallium atom.		Core = 1s, 2s2p, 3s3p	ECP guess small core   	7 UP, 6 DN   
2 2
0.716   5
4.165   2
0.712   5
4.325   1
Ge_ECP_SC      	32 	18	28	2.305	Germanium atom. 	Core = 1s, 2s2p, 3s3p	ECP guess small core	8 UP, 6 DN
2 2
0.670   5
3.334   3
0.676   5
3.296   1
As_ECP_SC      	33	18 	28	2.268	Arsen atom.		Core = 1s, 2s2p, 3s3p	ECP guess small core	9 UP, 6 DN
2 2
0.623   5
3.202   4
0.628   5
2.953   1
Se_ECP_SC      	34	18 	28	2.192	Selenium atom.		Core = 1s, 2s2p, 3s3p	ECP guess small core	9 UP, 7 DN
2 2 
0.582   5
3.192   4
0.590   5
2.910   2
Br_ECP_SC      	35	18 	28	2.154	Bromium atom.		Core = 1s, 2s2p, 3s3p	ECP guess small core   	9 UP, 8 DN
2 2
0.556   5
2.647   4
0.564   5
2.634   3
Kr_ECP_SC      	36	18 	28	2.116	Krypton atom.		Core = 1s, 2s2p, 3s3p	ECP guess small core   	9 UP, 9 DN    
2 2 
0.523   5
2.475   4
0.523   5
2.475   4
Krass      	36	18  	28	2.116	Krypton atom.        	Core = 1s,2s2p,3s3p     ECP guess small core    13 UP, 13 DN
2 2
0.523   7
2.475   6
0.523   7
2.475   6
    '''
    f = open(path+'xx_database_xx','w')
    f.write(data) 
    f.close() 
