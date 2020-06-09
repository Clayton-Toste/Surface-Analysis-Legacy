# ***************************************************************************
# Contains index, RGB triplets and a name for the color
# Any color index that is not set gets set to values for black -  0,0,0
# and the name is not set.
#
# Color values for R,G or B can be set to values from 0 up to 255
# Indices are from 1 to 128.  
#
# Lines with a "#" in the first column are ignored
# Lines are limited to 100 characters
#
# Herc Silverstein, Dec 1997
#    Feb 2000 	red ramp added 
#    		grad* ramp renamed blue*
# ***************************************************************************

# *** DO NOT change the names or indices for colors 1-96!
# *** You may change the rgb values for indices 1-96.
# *** You may change any of the fields for color 97 or greater.
# *** You may add, delete or change colors for index 97 or greater

#Index	RGB triplet 	Name
#-----	-----------	----

#  *********  These mimic the M1 color scheme ************
1	  0   0   0     black
2	160 160 160     gray
3	  0  0  180     dark blue
4	 30 30  225     blue
5	100 100 225     light blue
6	112 219 147     aquamarine
7	173 234 234     turquoise
8	  0 255 127     spring green
9	  0 100   0     dark green
10	 30 225  30     green
11	 50 204  50     lime green
12	153 204  30     yellow green
13	225 225  30     yellow
14	234 130  50     orange
15	142  35 107     maroon
16	225  30  30     red
17	255 152 163     pink
18	234 173 234     plum
#       purple is really magenta
19	225  30 225     purple
#	blue_purple is really blue violet
20	159  95 159     blue purple
21	255 255 255     white

# ************* M2 color additions in alphabetical order
22	165  42  42    brown
23	225 127  80    coral
24	105 105 105    dim gray
25	225 193  37    goldenrod
26	225 105 180    hot pink
27	107 142  35    olive
28	205 133  63    peru
29	160  82  45    sienna
30	 70 130 180    steel blue
31	216 191 216    thistle
32	245 222 179    wheat
			
# ************** Gradated colors.  
# ************** Blue ramp
33	  7   7 255    blue1
34	 15  15 255    blue2
35	 23  23 255    blue3
36	 31  31 255    blue4
37	 39  39 255    blue5
38	 47  47 255    blue6
39	 55  55 255    blue7
40       63  63 255    blue8
41       71  71 255    blue9
42       79  79 255    blue10
43       87  87 255    blue11
44       95  95 255    blue12
45      103 103 255    blue13
46      111 111 255    blue14
47      119 119 255    blue15
48      127 127 255    blue16
49      135 135 255    blue17
50      143 143 255    blue18
51      151 151 255    blue19
52      159 159 255    blue20
53      167 167 255    blue21
54      175 175 255    blue22
55      183 183 255    blue23
56      191 191 255    blue24
57      199 199 255    blue25
58      207 207 255    blue26
59      215 215 255    blue27
60      223 223 255    blue28
61      231 231 255    blue29
62      239 239 255    blue30
63      247 247 255    blue31
64      255 255 255    blue32

# ************** Red ramp
65   255   7   7  red1
66   255  15  15  red2
67   255  23  23  red3
68   255  31  31  red4
69   255  39  39  red5
70   255  47  47  red6
71   255  55  55  red7
72   255  63  63  red8
73   255  71  71  red9
74   255  79  79  red10
75   255  87  87  red11
76   255  95  95  red12
77   255 103 103  red13
78   255 111 111  red14
79   255 119 119  red15
80   255 127 127  red16
81   255 135 135  red17
82   255 143 143  red18
83   255 151 151  red19
84   255 159 159  red20
85   255 167 167  red21
86   255 175 175  red22
87   255 183 183  red23
88   255 191 191  red24
89   255 199 199  red25
90   255 207 207  red26
91   255 215 215  red27
92   255 223 223  red28
93   255 231 231  red29
94   255 239 239  red30
95   255 247 247  red31
96   255 255 255  red32

# ************** Some user defined colors.  Feel free to change these.
# *** dark grey-blue
97       44  60  60     user1
# *** navy
98        0   0 128     user2
# *** forestgreen
99       34 139  34     user3
# *** cyan
100       0 255 255     user4
# *** crimson 
101      220  20  60    user5
# *** chocolate
102      210 127  36    user6

# ************** Any other colors.  Don't go above 128.
103      192 192 192     user7
104      224 224 224     user8
105      255 30  0       user9
106      255 60  0       user10
107      255 90  0       user11
108      255 120 0       user12
109      255 150 0       user13
110      255 180 0       user14
111      255 210 0       user15
112      255 240 0       user16
113      255 255 0       user17
114      204 255 0       user18
115      153 255 0       user19
116      102 255 0       user20
117      51  255 0       user21
118      0   204 51      user22
119      0   153 102     user23
120      0   102 153     user24
121      0   51  204     user25
122      0   0   255     user26
123      15  0   255     user27
124      30  0   255     user28
125      45  0   255     user29
126      60  0   255     user30
127      75  0   255     user31
128      90  0   255     user32
129	255	0	31	userE
130	223	63	0	userD
131	191	127	25	userC
132	255	191	51	userT
133	223	223	0	userS
134	191	255	51	userG
135	159	255	0	userA
136	127	191	25	userP
137	95	255	51	userV
138	63	191	76	userI
139	31	255	102	userL
140	0	191	127	userM
141	159	255	191	userF
142	159	210	255	userW
143	255	223	159	userY
144	191	159	255	userH
145	95	63	255	userR
146	0	31	191	userK
147	191	0	191	userQ
148	255	0	159	userN
149	  0   255        51     user53
150	  0   255       102     user54
151	  0   255       153     user55
152	  0   255       204     user56
153	  0   204       255     user57
154	  0   153       255     user58
155	  0   102       255     user59
156	  0    51       255     user60
157	120     0       255     user61
158	180     0       255     user62
159	  0   0   38     user63
160	  0   0   0     user64
161	  0   0   0     user65
162	  0   0   0     user66
163	  0   0   0     user67
164	  0   0   0     user68
165	  0   0   0     user69
166	  0   0   0     user70
167	  0   0   0     user71
168	  0   0   0     user72
169	  0   0   0     user73
170	  0   0   0     user74
171	  0   0   0     user75
172	  0   0   0     user76
173	  0   0   0     user77
174	  0   0   0     user78
175	  0   0   0     user79
176	  0   0   0     user80
177	  0   0   0     user81
178	  0   0   0     user82
179	  0   0   0     user83
180	  0   0   0     user84
181	  0   0   0     user85
182	  0   0   0     user86
183	  0   0   0     user87
184	  0   0   0     user88
185	  0   0   0     user89
186	  0   0   0     user90
187	  0   0   0     user91
188	  0   0   0     user92
189	  0   0   0     user93
190	  0   0   0     user94
191	  0   0   0     user95
192	  0   0   0     user96
193	  0   0   0     user97
194	  0   0   0     user98
195	  0   0   0     user99
196	  0   0   0     user100
197	  0   0   0     user101
198	  0   0   0     user102
199	  0   0   0     user103
200	  0   0   0     user104
201	  0   0   0     user105
202	  0   0   0     user106
203	  0   0   0     user107
204	  0   0   0     user108
205	  0   0   0     user109
206	  0   0   0     user110
207	  0   0   0     user111
208	  0   0   0     user112
209	  0   0   0     user113
210	  0   0   0     user114
211	  0   0   0     user115
212	  0   0   0     user116
213	  0   0   0     user117
214	  0   0   0     user118
215	  0   0   0     user119
216	  0   0   0     user120
217	  0   0   0     user121
218	  0   0   0     user122
219	  0   0   0     user123
220	  0   0   0     user124
221	  0   0   0     user125
222	  0   0   0     user126
223	  0   0   0     user127
224	  0   0   0     user128
225	  0   0   0     user129
226	  0   0   0     user130
227	  0   0   0     user131
228	  0   0   0     user132
229	  0   0   0     user133
230	  0   0   0     user134
231	  0   0   0     user135
232	  0   0   0     user136
233	  0   0   0     user137
234	  0   0   0     user138
235	  0   0   0     user139
236	  0   0   0     user140
237	  0   0   0     user141
238	  0   0   0     user142
239	  0   0   0     user143
240	  0   0   0     user144
241	  0   0   0     user145
242	  0   0   0     user146
243	  0   0   0     user147
244	  0   0   0     user148
245	  0   0   0     user149
246	  0   0   0     user150
247	  0   0   0     user151
248	  0   0   0     user152
249	  0   0   0     user153
250	  0   0   0     user154
251	  0   0   0     user155
252	  0   0   0     user156
253	  0   0   0     user157
254	  0   0   0     user158
255	  0   0   0     user159
256	  0   0   0     user160

