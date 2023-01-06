import numpy as np
import pandas as pd
import random
from tqdm import tqdm # progress bar
import pickle # save dictionary in case there are some unlabeled courses


# Input:
# courses - list, all possible courseID
# days - list, [M, T, W, Th, F]
# profs - dict, key: prof, val: 2D list of time slots (time*days), if busy, entry = 1, else 0
# classrooms - list, all possible classrooms
# times - list, all possible time slots
# sems - 0: full, 1: first, 2: second
# acc - iterations
# M_init - previous Ms 
# (Step 1: first work on full semster courses -> 
# Step 2: first half based on M_init that we get from step 1 ->
# Step 3: second half based on M_init that we get from step 1)
# Therefore time slots & classrooms of full sems courses will be consistent
# 
# Output:
# schedule - [dataframe, dim(times * classrooms), courseID/break entries] for 5 days in a week
#
# greedy1: 
# choose the course with the most conflict only in the first iteration
# Since conflict matrix is very sparse, we assume conflicted courses are assigned
# in the FIRST iteration of "for d in tqdm(range(len(days)), desc=" days", position=0):"
# and don't follow the priority from the second iteration through shuffling.
def schedule(courses, days, profs, classrooms, times, sems, acc=0, M_init=None):
    if M_init != None:
        M_all = M_init
    else:
        M_all = [pd.DataFrame(np.zeros((len(times), len(classrooms)), dtype=np.uint8), 
                            columns=[classrooms], index=[times]) for i in range(len(days))] # dataframes for Mon-Fri
    for d in tqdm(range(len(days)), desc=" days", position=0):
        initr, initc = (0,0)
        labeled = []
        unlabeled = []
        for course in tqdm(courses, desc=" courses", position=1, leave=False): 
            # ignore pointr, pointc. No use.
            new_pointr, new_pointc, new_M_all, new_prof_dict, islabeled = assign(course, days[d], profs, M_all, initr, initc)
            if not islabeled:
                unlabeled.append(course)
            else:
                labeled.append(course)
            initr, initc, M_all, profs = new_pointr, new_pointc, new_M_all, new_prof_dict
        courses = set(unlabeled.copy())
        if acc > 0:
            random.shuffle(list(courses))
    for m_idx in range(len(M_all)):
        M_all[m_idx].to_csv(f'S22Schedule_sems{sems}_{m_idx}.csv')
    if acc >= 5: # Don't really enter this branch 
        print("Too many recursion. Find a better solution.")
        return profs, M_all, unlabeled
    # reassign
    elif len(unlabeled) > 0: # Don't really enter this branch
        unlabled_cp = list(set(unlabeled[:]))
        random.shuffle(unlabled_cp)
        return schedule(unlabled_cp, days, profs, classrooms, times, sems, acc+1, M_all)
    else: # Base Case, len(unlabeled) == 0
        print("Scheduling problem solved.")
    return profs, M_all, None


# constraints (See helper functions for details)
# 1) capacity 
# 2) prof cannot teach two classes at the same time slot
# 3) one major's required courses should not be at the same time slot in any classroom
#
# if class also on another day assign to same time on that day
#
# Input:
# course - Lecture type course
# day - a day from Mon to Fri
# prof_dict - prof schedules
# M_all - 5 dataframes for schedules
# pointr, pointc - ignore. No use.
# 
# Output: 
# M_all - updated schedule
# prof_dict - updated prof schedule
# last boolean arg: True/False, if the course is assigned -> True, o.w. -> False.
#
# greedy2: 
# choose the classroom with the lowest capacity - sort classroom by capacity when initializing M, loop through classrooms
# choose the earliest timeslot - loop through time slot from the beginning  
def assign(course, day, prof_dict, M_all, pointr, pointc):
    crseNum = course.crseNum
    sect = course.sect
    prof = course.prof
    duration = course.duration
    days = course.days


    initr = 0
    initc = 0
    rows, cols = M_all[day].shape
    roomidx = initc
    rooms = list(M_all[day].columns)
    timeidx = initr
    
    if course.prof is not None: # Ignore TAs and TBA professors
        profs = course.prof.split(",") # if there are multiple professors, assume all profs teach this course
    else:
        profs = []

    while initc <= roomidx < cols: # loop through all rooms from small classroom to big ones, rooms already sorted
        while initr <= timeidx < rows: # loop though all times starting from 8am.
            times = range(timeidx, timeidx+(int(duration)//5)) 
            
            # constraints
            if timeidx+(int(duration)//5) < rows and notFilled(times, rooms, roomidx, days, M_all) \
                and capacity(course,rooms[roomidx]) \
                and schedule_prof(times,day,profs,prof_dict) \
                        and major(course,times,M_all):
                # all valid days
                valid_days = days 
                days_dict = {'M':0,'T':1,'W':2,'R':3,'F':4}               
                for d in valid_days:
                    for prof in profs:
                        prof_dict[prof][times,days_dict[d]] = [1]*len(times) # update prof schedules
                    M_all[days_dict[d]].loc[times,rooms[roomidx]] = [crseNum + '-' + sect]*len(times)      
                    # add 15 minutes break after a course
                    if timeidx+(int(duration)//5) + 3 > rows:
                        M_all[days_dict[d]].loc[range((timeidx+(int(duration)//5)),rows),rooms[roomidx]] = ['break']*(rows-(timeidx+(int(duration)//5)))
                    else:
                        M_all[days_dict[d]].loc[range(timeidx+(int(duration)//5),timeidx+(int(duration)//5)+3),rooms[roomidx]] = ['break']*3   
                pointr = timeidx # ignore
                pointc = roomidx # ignore
                # for debugging
                #print("assigning ",crseNum+'-'+sect, 'to',rooms[roomidx],'at',times)
                return pointr, pointc, M_all, prof_dict, True
            timeidx += 1
        timeidx = pointr # ignore
        roomidx += 1
    return initr, initc, M_all, prof_dict, False

#conflict checks helpers

#ensures enough capacity, max enrollment <= seats 
def capacity(course, classroom): 
    return (int(course.maxEnr) \
            <= classrooms[classrooms['Rooms'].str.contains(classroom[0])].loc[:,'SEATS'].values[0])

#prof not teaching another class at the same time
def schedule_prof(times, day, profs, profs_dict): 
    for prof in profs:
        for t in times:
            t_max,_ = profs_dict[prof].shape
            if t >= t_max or profs_dict[prof][t,day] == 1:
                return False
    return True

#another class in some major's pathway not at same time on any other day
def major(course, times, M_all): 
    conflicts = pd.DataFrame(major_conflicts[course.crseNum]) #get correct row
    conflicts = conflicts.loc[(conflicts[course.crseNum] != 0),:] #only take cols where conflict
    conflicts = set(conflicts.columns) #column names are courses with conflict
    prev_scheduled = []
    for M_d in M_all:
        prev_scheduled.extend(list(M_d.loc[times,:].values.flatten().tolist()))
    return len(set(prev_scheduled).intersection(conflicts)) == 0
        
# make sure all candidate time slots are empty
def notFilled(time_interval, rooms, roomidx, days, M_all):
    days_dict = {'M':0,'T':1,'W':2,'R':3,'F':4}  
    for d in days:
        tar = np.array(M_all[days_dict[d]].loc[time_interval,rooms[roomidx]]) 
        if np.any(tar != 0):
            return False
    return True

# all global variables here
if __name__ == "__main__":
    #random.seed(21393) # change this number to get different results in shuffle
    
    major_conflicts = pd.read_excel('filled_conflicts_undergrad_fin.xlsx') # conflict matrix for undergrad courses
    major_conflicts.dropna(axis=0,how='any',inplace=True)
    major_conflicts.columns = major_conflicts.columns.astype(str)

    classrooms = pd.read_excel('S22_Registrar_Classrooms_1nov21.xlsx') # all classrooms, with number of seats
    S22 = pd.read_excel('S22_Registrar_Schedule_Courses_1nov21.xlsx') # provided spring sems schedule
    
    # data preprocessing, convert everything to string
    full = pd.ExcelFile("output_full.xlsx") 
    monF = pd.read_excel(full, 'Full_Mon', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    tueF = pd.read_excel(full, 'Full_Tue', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    wedF = pd.read_excel(full, 'Full_Wed', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    thurF = pd.read_excel(full, 'Full_Thur', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    friF = pd.read_excel(full, 'Full_Fri', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    mini = pd.ExcelFile("output_mini.xlsx")
    monFH = pd.read_excel(mini, 'fh_Mon', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    tueFH = pd.read_excel(mini, 'fh_Tue', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    wedFH = pd.read_excel(mini, 'fh_Wed', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    thurFH = pd.read_excel(mini, 'fh_Thur', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    friFH = pd.read_excel(mini, 'fh_Fri', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})

    monSH = pd.read_excel(mini, 'sh_Mon', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    tueSH = pd.read_excel(mini, 'sh_Tue', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    wedSH = pd.read_excel(mini, 'sh_Wed', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    thurSH = pd.read_excel(mini, 'sh_Thur', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    friSH = pd.read_excel(mini, 'sh_Fri', converters={'CRSE#': lambda x: str(x), 
                                                      'MAX ENR': lambda x: str(x), 
                                                      'Duration': lambda x: str(x)})
    

# define Lecture type
class Lecture: 
    def __init__(self, crseNum, sect, maxEnr, prof, duration, days): 
        self.crseNum = crseNum.split("/")[0]
        self.sect = sect
        self.maxEnr = maxEnr
        self.prof = prof
        self.duration = duration
        self.days = days

    def __eq__(self,other): # make Lecture hashable, since set operations are needed later
        return (self.crseNum == other.crseNum) and \
        (self.sect == other.sect) 

    def __hash__(self):  # make Lecture hashable, since set operations are needed later
        return hash(self.crseNum + self.sect)

def listLectures(classDay):
    lst = np.array(classDay)
    output = [] 
    for i in range(len(lst)):
        crseNum = lst[i][3]
        sect = lst[i][4]
        maxEnr = lst[i][9]
        if (pd.isna(lst[i][14])):
            prof = None
        else:
            prof = lst[i][14]
        duration = lst[i][21]  
        day = lst[i][11]    
        output.append(Lecture(crseNum, sect, maxEnr, prof, duration, day))
    return output

mF =  listLectures(monF)
tF =  listLectures(tueF)
wF =  listLectures(wedF)
rF =  listLectures(thurF)
fF =  listLectures(friF)

mFH = listLectures(monFH)
tFH = listLectures(tueFH)
wFH = listLectures(wedFH)
rFH = listLectures(thurFH)
fFH = listLectures(friFH)

mSH = listLectures(monSH)
tSH = listLectures(tueSH)
wSH = listLectures(wedSH)
rSH = listLectures(thurSH)
fSH = listLectures(friSH)

full_schedules = [mF, tF, wF, rF, fF]
first_schedules = [mFH, tFH, wFH, rFH, fFH]
second_schedules = [mSH, tSH, wSH, rSH, fSH]

# contain all Lecture type courses for full/first/second sems
# So that we can get some attributes easily later
full_courses = []
first_courses = []
second_courses = []

for i in range(len(full_schedules)):
  for j in range(len(full_schedules[i])):
    full_courses.append(full_schedules[i][j])

for i in range(len(first_schedules)):
  for j in range(len(first_schedules[i])):
    first_courses.append(first_schedules[i][j])

for i in range(len(second_schedules)):
  for j in range(len(second_schedules[i])):
    second_courses.append(second_schedules[i][j])

# all undergrad or undergrad/grad crseNum (courseID), exclude courses only for graduates
courseID = list(filter(lambda x: x.isdigit(), set(np.array(major_conflicts.columns, dtype=str)))) 

full_courseID = []
for cfu in full_courses:
    if cfu.crseNum in courseID:
        full_courseID.append(cfu.crseNum)

first_courseID = []
for cf in first_courses:
    if cf.crseNum in courseID:
        first_courseID.append(cf.crseNum)

second_courseID = []
for cs in second_courses:
    if cs.crseNum in courseID:
        second_courseID.append(cs.crseNum)

# sort full/first half/second half courses by priorities
full_priority = list(major_conflicts.loc[:,full_courseID].sum(axis=0)) # sum of all conflicts with this course
full_courses_with_priority = sorted(zip(full_courses,full_priority),key=lambda x: x[1], reverse=True)
full_courses = list(map(lambda x: x[0],full_courses_with_priority))

first_priority = list(major_conflicts.loc[:,first_courseID].sum(axis=0)) # sum of all conflicts with this course
first_courses_with_priority = sorted(zip(first_courses,first_priority),key=lambda x: x[1], reverse=True)
first_courses = list(map(lambda x: x[0],first_courses_with_priority))

second_priority = list(major_conflicts.loc[:,second_courseID].sum(axis=0)) # sum of all conflicts with this course
second_courses_with_priority = sorted(zip(second_courses,second_priority),key=lambda x: x[1], reverse=True)
second_courses = list(map(lambda x: x[0],second_courses_with_priority))

days = range(5) # [M, T, W, TR, F] 

hours = 13 #total allowable hours for scheduling
times = [i for i in range(hours * 60 // 5)] #t in [0, n] represent evenly spaced blocks of 5 min

# set up prof schedules as dictionaries
prof_names = list(set(np.array(S22['INSTRUCTOR(S)'])))
profs = dict(zip(prof_names, [ np.zeros((len(times), len(days)),dtype=np.uint8) for i in range(len(prof_names))]))

# sort classrooms by capacity
classrooms["Rooms"] = classrooms.loc[:, 'BLDG_SIS'] + classrooms.loc[:, 'SISROOMS_NUMBER'] # buildings + room number
room_with_capacity = sorted(zip(list(classrooms["Rooms"]),list(classrooms['SEATS'])),key=lambda x: x[1])
classrooms_L = list(map(lambda x: x[0],room_with_capacity))

# full sems
print("Working on the full semester...")
profs_full, M_full, unlabeled_full = schedule(set(full_courses), days, profs, classrooms_L, times, 0, 0)

# first sems
print("Working on the first half of semester...")
_,_, unlabeled_first = schedule(set(first_courses), days, profs_full, classrooms_L, times, 1, 0, M_init = M_full)

# second sems
print("Working on the second half of semester...")
_,_, unlabeled_second = schedule(set(second_courses), days, profs_full, classrooms_L, times, 2, 0, M_init = M_full)

# save unlabeled courses
# All of them are empty finally
with open('unlabeled_full.pickle', 'wb') as handle:
    pickle.dump(unlabeled_full, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('unlabeled_first.pickle', 'wb') as handle:
    pickle.dump(unlabeled_first, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('unlabeled_second.pickle', 'wb') as handle:
    pickle.dump(unlabeled_second, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Warning: Don't open csv/excel files while running your code!!!!!!!!
