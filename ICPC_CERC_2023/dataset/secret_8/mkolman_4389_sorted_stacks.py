from dataclasses import dataclass
class Stack:
    def __init__(self, idx, s):
        self.idx = idx
        self.c, *self.s = s
        while self.s and self.s[-1] == 0:
            self.s.pop()
    
    def take(self):
        return self.s.pop()

    def put(self, x):
        self.s.append(x)
        if len(self.s) > self.c:
            raise ValueError("Stack is full")
    
    def empty(self):
        return self.c - len(self.s)

def find(stacks, x):
    for i, stack in enumerate(stacks):
        if x in stack.s:
            return i, stack.s.index(x)
    raise ValueError("Value not found")

def swap(stacks, src_idx, dst_idx):
    if src_idx == dst_idx:
        return
    stacks[dst_idx].put(stacks[src_idx].take())
    print(src_idx+1, dst_idx+1)

def clear(stacks, stack_idx, num, from_stack):
    """Take `num` cards from `stack_idx` and put them on any other stack `from_stack` or higher."""
    if num == 0:
        return 0
    for to_stack in range(from_stack, len(stacks)):
        if to_stack == stack_idx:
            continue
        while stacks[to_stack].empty() > 0 and num > 0:
            swap(stacks, stack_idx, to_stack)
            num -= 1
        if num == 0:
            return 0
    return num

def move(stacks, src, dst):
    """Move src to dst without disturbing the order of stacks dst and before."""
    if src == dst:
        return
    # Move as much junk from src to garbage stacks as possible
    top = len(stacks[src[0]].s) - src[1] - 1 # number of cards above src
    top = clear(stacks, src[0], top, dst[0]+1)
    if top:
        # Now all stacks other than dst and src are full.
        # Choose any stack that is not src or dst as tmp
        tmp = len(stacks) - 1
        while tmp in [src[0], dst[0]]:
            tmp -= 1
        # Clear the top of tmp to use for moving src
        swap(stacks, tmp, dst[0])
        # Move the top of src to tmp
        for _ in range(top):
            swap(stacks, src[0], dst[0])
        # Move src to tmp
        swap(stacks, src[0], tmp)
        # Move values above src back to their original place
        for _ in range(top):
            swap(stacks, dst[0], src[0])
        # Move src back to the TOP of its original place
        swap(stacks, tmp, src[0])
        # Restore tmp to its original state
        swap(stacks, dst[0], tmp)
    # Move src to dst
    swap(stacks, src[0], dst[0])

def main():
    n, k = map(int, input().split())
    stacks = []
    for i in range(k):
        stacks.append(Stack(i+1, map(int, input().split())))

    # Fill the stacks in order
    i = 0
    for stack_idx in range(k):
        clear(stacks, stack_idx, len(stacks[stack_idx].s), stack_idx+1)
        for idx in range(stacks[stack_idx].c):
            i += 1
            if i > n:
                print(0, 0)
                return
            move(stacks, find(stacks, i), (stack_idx, idx))

if __name__ == '__main__':
    main()