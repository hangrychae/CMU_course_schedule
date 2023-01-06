import pandas as pd
import numpy as np

S22 = pd.read_excel('final/S22_Registrar_Schedule_Courses_1nov21.xlsx')
S22_courses = np.unique(np.array(S22.loc[(S22["Course Lvl"] == "U")|(S22["Course Lvl"] == "U/G"),"CRSE#"]))
S22_courses = list(map(lambda x: x.split("/")[0], S22_courses))
df = pd.read_excel("final/Major_Conflicts.xlsx",sheet_name = 'Sheet1',converters={'CRSE#': lambda x: ('0' + str(x) if len(str(x)) == 4 else str(x))})
del df["CLG"]
del df["DEPT"]
del df["TITLE"]
df.set_index('CRSE#',inplace=True)
df = df.fillna(0)
df_p = pd.read_csv("undergrad_priority_output.csv")
df_p.dropna(how='any',inplace=True)
rows,cols = df_p.shape
print(df.shape)
#import pdb; pdb.set_trace()
for r in range(rows):
    if df.index[r] in S22_courses:
        edge = df_p.iloc[r,0]
        priority = df_p.iloc[r,1]
        A = edge.split(",")[0][1:]
        B = edge.split(",")[1][:-1]
        A = A[1:3] + A[4:-1]
        B = B[2:4] + B[5:-1]
        df.loc[A,B] = priority
        df.loc[B,A] = priority
df.dropna(how='all',inplace=True)
df.to_excel("final/filled_conflicts_undergrad_fin.xlsx")