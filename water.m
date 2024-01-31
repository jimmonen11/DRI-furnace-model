function [rho_H2O, mu_H2O, k_H2O, C_H2O] = water(T_H2O)
%CoolProp water properties at 1.15 bara as function of T

properties = [...
290	998.810957846585	0.00108396587239529	    0.592307321025754	4186.54962561573
300	996.563758825863	0.000853741110538816	0.609508385411362	4180.59281577039
310	993.390323903581	0.000693330518053731	0.624277902175849	4179.20302393507
320	989.433474507222	0.000576729041693834	0.637003698122406	4180.49908988052
330	984.793411795782	0.000489151035887246	0.647919146202427	4183.61758994246
340	979.542760740358	0.000421637438943552	0.657175897646496	4188.25993979176
350	973.735223433550	0.000368473827067594	0.664882165497081	4194.43410345288
360	967.410814235632	0.000325859738863381	0.671122912613493	4202.30441988547
370	960.599140146311	0.000291179175134149	0.675970297803990	4212.10793043996
380	0.675438261236059	1.24820639994272e-05	0.0252028085426623	2082.96839985096
390	0.656779857886044	1.28721780484354e-05	0.0260416757204778	2049.39241175825
400	0.639284130523932	1.32647942148055e-05	0.0268914135756934	2027.15100173421
410	0.622812870293532	1.36597537524474e-05	0.0277528358417335	2011.40514334885
420	0.607257641418547	1.40568962176074e-05	0.0286264958591902	2000.14101041670
430	0.592529213070273	1.44560630655696e-05	0.0295127444171416	1992.21908905625
440	0.578552391982740	1.48570995428409e-05	0.0304117775506470	1986.87364857678
450	0.565262808161383	1.52598557811270e-05	0.0313236737052506	1983.54988617955
460	0.552604678160489	1.56641874270011e-05	0.0322484222912654	1981.83220837680
470	0.540529166888444	1.60699559756726e-05	0.0331859456008890	1981.40368309932
480	0.528993152728450	1.64770289095973e-05	0.0341361156564940	1982.01964002676
490	0.517958277205058	1.68852797068726e-05	0.0350987671867471	1983.48930567318
500	0.507390200878015	1.72945877627547e-05	0.0360737076366436	1985.66261595817
510	0.497258011693236	1.77048382537711e-05	0.0370607248980988	1988.42057456395
520	0.487533747957169	1.81159219647244e-05	0.0380595932830753	1991.66811399791
530	0.478192008818588	1.85277350926853e-05	0.0390700781375361	1995.32875436200
540	0.469209632513098	1.89401790378202e-05	0.0400919394012950	1999.34056948687
550	0.460565427782696	1.93531601879482e-05	0.0411249343483119	2003.65311348344
560	0.452239947544570	1.97665897016494e-05	0.0421688196884036	2008.22505927541
570	0.444215296519126	2.01803832932812e-05	0.0432233531704848	2013.02236948529
580	0.436474966449078	2.05944610222102e-05	0.0442882947961643	2018.01686863715
590	0.429003693960088	2.10087470878193e-05	0.0453634077284798	2023.18512029433
600	0.421787337173110	2.14231696313061e-05	0.0464484589620027	2028.50753769549
610	0.414812767979058	2.18376605449007e-05	0.0475432198061862	2033.96767454956
620	0.408067777497815	2.22521552888470e-05	0.0486474662226773	2039.55165587281
630	0.401540992715398	2.26665927162896e-05	0.0497609790486220	2045.24771848010
640	0.395221802660839	2.30809149060618e-05	0.0508835441311999	2051.04583795037
650	0.389100292773802	2.34950670032684e-05	0.0520149523933031	2056.93742426238
660	0.383167186343716	2.39089970674832e-05	0.0531549998460858	2062.91507233114
670	0.377413792085319	2.43226559283328e-05	0.0543034875608199	2068.97235672378
680	0.371831957064141	2.47359970482069e-05	0.0554602216098854	2075.10366214954
690	0.366414024306424	2.51489763918130e-05	0.0566250129846658	2081.30404308980
700	0.361152794527113	2.55615523022860e-05	0.0577976774964877	2087.56910729714
710	0.356041491491406	2.59736853835576e-05	0.0589780356654444	2093.89491894762
720	0.351073730593311	2.63853383886922e-05	0.0601659126009164	2100.27791805206
730	0.346243490291469	2.67964761139002e-05	0.0613611378767789	2106.71485337677
740	0.341545086090235	2.72070652979481e-05	0.0625635454036332	2113.20272663179
750	0.336973146794350	2.76170745266921e-05	0.0637729732998802	2119.73874608850
760	0.332522592799784	2.80264741424733e-05	0.0649892637630358	2126.32028811015
770	0.328188616212544	2.84352361581220e-05	0.0662122629423609	2132.94486533817
780	0.323966662612265	2.88433341753316e-05	0.0674418208136125	2139.61010048614
790	0.319852414298970	2.92507433071714e-05	0.0686777910565150	2146.31370486419
800	0.315841774879944	2.96574401045203e-05	0.0699200309353823	2153.05346089584
810	0.311930855069817	3.00634024862146e-05	0.0711684011831922	2159.82720800442
820	0.308115959590965	3.04686096727124e-05	0.0724227658893092	2166.63283134167
830	0.304393575073578	3.08730421230880e-05	0.0736829923909690	2173.46825291013
840	0.300760358865436	3.12766814751806e-05	0.0749489511685769	2180.33142469787
850	0.297213128670854	3.16795104887286e-05	0.0762205157448209	2187.22032349982
860	0.293748852946510	3.20815129913336e-05	0.0774975625875605	2194.13294714741
870	0.290364641989167	3.24826738271027e-05	0.0787799710164273	2201.06731190847
880	0.287057739656774	3.28829788078285e-05	0.0800676231130443	2208.02145085358
890	0.283825515670119	3.32824146665741e-05	0.0813604036347604	2214.99341301441
900	0.280665458447314	3.36809690135357e-05	0.0826581999317726	2221.98126318486
910	0.277575168427906	3.40786302940649e-05	0.0839609018674969	2228.98308223732
920	0.274552351847449	3.44753877487366e-05	0.0852684017420170	2235.99696784507
930	0.271594814926970	3.48712313753582e-05	0.0865805942183688	2243.02103551821
940	0.268700458444990	3.52661518928176e-05	0.0878973762512203	2250.05341987406
950	0.265867272662681	3.56601407066764e-05	0.0892186470167080	2257.09227607546
960	0.263093332575299	3.60531898764175e-05	0.0905443078375716	2264.13578138052
970	0.260376793465419	3.64452920842642e-05	0.0918742619995668	2271.18213675676
980	0.257715886735564	3.68364406054877e-05	0.0932084155251698	2278.22956852024
990	0.255108915999774	3.72266292801298e-05	0.0945466756872654	2285.27632996697
1000	0.252554253415314	3.76158524860688e-05	0.0958889517310973	2292.32070297031
1010	0.250050336237357	3.80041051133585e-05	0.0972351549208067	2299.36099952239
1020	0.247595663580818	3.83913825397799e-05	0.0985851983419077	2306.39556320263
1030	0.245188793374840	3.87776806075414e-05	0.0999389968519848	2313.42277055987
1040	0.242828339496570	3.91629956010731e-05	0.101296467033223	2320.44103239799
1050	0.240512969071950	3.95473242258583e-05	0.102657527146686	2327.44879495750
1060	0.238241399932169	3.99306635882527e-05	0.104022097088255	2334.44454098828
1070	0.236012398215351	4.03130111762420e-05	0.105390098346163	2341.42679071024
1080	0.233824776103805	4.06943648410910e-05	0.106761453960033	2348.39410266068
1090	0.231677389687935	4.10747227798423e-05	0.108136088481375	2355.34507442843
1100	0.229569136948563	4.14540835186207e-05	0.109513927935447	2362.27834327592
1110	0.227498955850011	4.18324458967073e-05	0.110894899784441	2369.19258665153
1120	0.225465822536901	4.22098090513426e-05	0.112278932891923	2376.08652259513
1130	0.223468749628077	4.25861724032266e-05	0.113665957488476	2382.95891004051
1140	0.221506784601603	4.29615356426796e-05	0.115055905138483	2389.80854901877
1150	0.219579008265155	4.33358987164350e-05	0.116448708708016	2396.63428076722
1160	0.217684533306578	4.37092618150308e-05	0.117844302333762	2403.43498774856
1170	0.215822502919710	4.40816253607742e-05	0.119242621392968	2410.20959358524
1180	0.213992089500945	4.44529899962499e-05	0.120643602474328	2416.95706291425
1190	0.212192493412294	4.48233565733482e-05	0.122047183349804	2423.67640116719
1200	0.210422941807000	4.51927261427872e-05	0.123453302947321	2430.36665428106
1210	0.208682687514040	4.55610999441064e-05	0.124861901324304	2437.02690834465
1220	0.206971007978079	4.59284793961103e-05	0.126272919642035	2443.65628918551
1230	0.205287204251669	4.62948660877394e-05	0.127686300140772	2450.25396190262
1240	0.203630600036713	4.66602617693509e-05	0.129101986115621	2456.81913034919
1250	0.202000540772383	4.70246683443879e-05	0.130519921893128	2463.35103657044
1260	0.200396392766890	4.73880878614219e-05	0.131940052808543	2469.84896020067
1270	0.198817542370647	4.77505225065481e-05	0.133362325183763	2476.31221782387
1280	0.197263395188542	4.81119745961206e-05	0.134786686305890	2482.74016230195
1290	0.195733375329161	4.84724465698096e-05	0.136213084406418	2489.13218207441
1300	0.194226924688964	4.88319409839675e-05	0.137641468640992	2495.48770043310
1310	0.192743502269507	4.91904605052886e-05	0.139071789069742	2501.80617477568
1320	0.191282583525954	4.95480079047503e-05	0.140503996638153	2508.08709584075
1330	0.189843659745202	4.99045860518225e-05	0.141938043158469	2514.32998692802
1340	0.188426237452063	5.02601979089326e-05	0.143373881291593	2520.53440310617
1350	0.187029837842030	5.06148465261759e-05	0.144811464529481	2526.69993041127
1360	0.185653996239252	5.09685350362597e-05	0.146250747177998	2532.82618503811
1370	0.184298261578400	5.13212666496704e-05	0.147691684340236	2538.91281252680
1380	0.182962195909227	5.16730446500540e-05	0.149134231900264	2544.95948694694
1390	0.181645373922641	5.20238723898009e-05	0.150578346507300	2550.96591008116
1400	0.180347382497228	5.23737532858254e-05	0.152023985560291	2556.93181060993
];

tempinterp = properties(:,1);
interp = interp1(tempinterp, properties, T_H2O);

rho_H2O = interp(2); %density, kg/m^3
mu_H2O = interp(3); %viscocity, Pa*s
k_H2O = interp(4); %conductive heat trans coef (W/mK)
C_H2O = interp(5); %specific heat capacity (J/kgK)
