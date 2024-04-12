print("happy birthday")

for number in range(0,50):
    print("{}: happy birthday".format(number))

# Challenge [Loops-1]: You are a robot that says the following 3 sentences over and over. Please write the transcript of your conversation.
sentences = ["I like your stoyl", "I think its floy", "I would never block yours"]

for number in range(0,12):
    print(sentences) 


# problem solved 
for number in range(0,12):
    for rizz in sentences:
        print(rizz)



