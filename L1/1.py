s = input("enter a string: ")

unique_letters = set()
for ch in s:
    if ch.isalpha():
        unique_letters.add(ch.upper())

# sort so output order is alphabetical
alphabet = "".join(sorted(unique_letters))

print("alphabet:", alphabet)
