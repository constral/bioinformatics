s = input("enter a string: ")

counts = [0] * 26
total = 0

for ch in s:
    cl = ch.lower()
    if 'a' <= cl <= 'z':
        idx = ord(cl) - ord('a')
        counts[idx] += 1
        total += 1

print("letter frequency:")
for i in range(26):
    if counts[i] > 0:
        pct = counts[i] * 100.0 / total
        print(f"{chr(ord('a') + i)} : {pct:.2f}%")
