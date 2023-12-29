n, k = map(int, input().split())
capacity, stack = [], []
pos, tar = {}, {}
prev = 0
for i in range(k):
    m, *st = map(int, input().split())
    while st and st[-1]==0: st.pop()
    for j,x in enumerate(st): pos[x] = (i,j)
    for j in range(m): tar[1+prev+j] = (i,j)
    prev += m
    stack.append(st)
    capacity.append(m)

# make a move, update positions
def move(i,j):
    print(1+i,1+j)
    x = stack[i].pop()
    pos[x] = (j, len(stack[j]))
    stack[j].append(x)

# move from src without filling stage to capacity and completely avoiding blocked
def remove(src, stage=None, blocked=None):
    for i in range(k):
        if i==src or i==blocked: continue
        if len(stack[i]) < capacity[i] - int(i==stage):
            move(src,i)
            return True
    return False

for it in range(2):  # keep one space in each stack in round 1
    for x in range(1,n+1):
        if x not in pos: continue
        i1, j1 = pos[x]
        i2, _ = tar[x]
        if it==0 and capacity[i2]==1: continue  # keep one space
        # make them different
        if i1 == i2:
            while len(stack[i1])-1 >= j1:  # disperse source stack
                remove(i1)
            i1, j1 = pos[x]
        # find stage stack
        for i3 in range(k):
            if capacity[i3] > 0 and i3 not in [i1, i2]: break
        else: i3 = i2  # special case with only two items left
        # clear stage stack
        if len(stack[i3]) == capacity[i3]:
            remove(i3,i3,i1)
        # disperse source stack above x without filling stage stack
        while len(stack[i1])-1 > j1:
            remove(i1,i3)
        # disperse target stack (avoid i1 where x is waiting)
        blocked = i1
        while len(stack[i2]) > 0:
            if not remove(i2,i3,blocked):  # unless forced to fill stage stack with x and free i1
                move(i1,i3)
                blocked = None
        # move x into position
        if blocked is None: move(i3,i2)  # x at top of stage stack
        else: move(i1,i2)  # x at top of source stack
        # reduce problem size
        stack[i2] = []
        capacity[i2] -= 1
        pos.pop(x)
print(0,0)
