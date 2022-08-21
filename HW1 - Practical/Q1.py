import string
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(1050)

class Huffman_Node:
    def __init__(self, key, frequency, right=None, left=None, code=''):
        self.key = key
        self.frequency = frequency
        self.right = right
        self.left = left
        self.code = code
        
def encoding_huffman_tree(inputs, flag=True):
    if flag:
        huffman_nodes = []
        for key, frequency in inputs.items():
            huffman_nodes.append(Huffman_Node(key, frequency))
        length = len(huffman_nodes)
        for i in range(length - 1):
            min_index = 0
            for j in range(1, length - i):
                if huffman_nodes[j].frequency < huffman_nodes[min_index].frequency:
                    min_index = j
            temp = huffman_nodes[length - i - 1]
            huffman_nodes[length - i - 1] = huffman_nodes[min_index]
            huffman_nodes[min_index] = temp
        inputs = huffman_nodes
    length = len(inputs)
    if length <= 1:
        return inputs[0]
    merged_huffman_node = Huffman_Node('merged', inputs[length - 1].frequency + inputs[length - 2].frequency, inputs[length - 1], inputs[length - 2], '')
    inputs.pop()
    inputs.pop()
    length = len(inputs)
    if length <= 0:
        inputs.append(merged_huffman_node)
    else:
        temp = length
        for i in range(length):
            if merged_huffman_node.frequency > inputs[i].frequency:
                temp = i
                break
        inputs.insert(i, merged_huffman_node)
    return encoding_huffman_tree(inputs, False)

def decoding_huffman_tree(huffman_node, alphabet_code_dict, flag=True):
    if huffman_node.right == None:
        if flag:
            alphabet_code_dict[huffman_node.key] = '0'
        else:
            alphabet_code_dict[huffman_node.key] = huffman_node.code
    else:
        huffman_node.right.code = huffman_node.code + '0'
        decoding_huffman_tree(huffman_node.right, alphabet_code_dict, False)
        huffman_node.left.code = huffman_node.code + '1'
        decoding_huffman_tree(huffman_node.left, alphabet_code_dict, False)

def huffman_compressor(sentence, alphabet_code_dict):
    sentence_binary_huffman = ''
    for i in sentence:
        sentence_binary_huffman += alphabet_code_dict[i]
    return sentence_binary_huffman

def n_p_sequences(n, p):
    previous_n_sequences_dict = {}
    for i in range(1, n + 1):
        if i == 1:
            previous_n_sequences_dict['0'] = 1 - p
            previous_n_sequences_dict['1'] = p
        else:
            next_n_sequences_dict = {}
            for i in previous_n_sequences_dict:
                next_n_sequences_dict[i + '0'] = previous_n_sequences_dict[i] * (1 - p)
                next_n_sequences_dict[i + '1'] = previous_n_sequences_dict[i] * p
            previous_n_sequences_dict = next_n_sequences_dict
    return previous_n_sequences_dict

def n_huffman_compressor(random_numbers_binary, n, all_n_p_sequences_code):
    random_numbers_binary_huffman = ''
    for i in range(int(len(random_numbers_binary) / n)):
        random_numbers_binary_huffman += all_n_p_sequences_code[random_numbers_binary[n * i: n * (i + 1)]]
    return random_numbers_binary_huffman

sentence = 'go go godzila'
sentence_ascii = [ord(char) for char in sentence]
sentence_binary = ''
for char_ascii in sentence_ascii:
    sentence_binary += bin(char_ascii)[2:].zfill(7)
print(sentence_ascii)
print(sentence_binary)
alphabet_string = string.ascii_lowercase + ' '
alphabet_dict = {}
for alphabet_char in alphabet_string:
    alphabet_dict[alphabet_char] = 0.0
n = len(sentence)
for char in sentence:
    alphabet_dict[char] += 1 / n
print(alphabet_dict)
huffman_tree_root = encoding_huffman_tree(alphabet_dict)
alphabet_code_dict = alphabet_dict
decoding_huffman_tree(huffman_tree_root, alphabet_code_dict)
print(alphabet_code_dict)
sentence_binary_huffman = huffman_compressor(sentence, alphabet_code_dict)
print(sentence_binary_huffman)
block_lengths = np.array(list(range(1, 11)))
repeat_number = (np.ceil(1000 / block_lengths) * block_lengths).astype(np.int)
lengths = np.zeros((3, 10), dtype=np.int)
index = 0
for p in [0.5, 0.75, 0.97]:
    random_numbers = np.random.rand(np.max(repeat_number))
    random_numbers_binary = ''
    for i in random_numbers:
        if i <= p:
            random_numbers_binary += '1'
        else:
            random_numbers_binary += '0'
    for n in block_lengths:
        all_n_p_sequences = n_p_sequences(n, p)
        huffman_tree_root = encoding_huffman_tree(all_n_p_sequences)
        all_n_p_sequences_code = all_n_p_sequences
        decoding_huffman_tree(huffman_tree_root, all_n_p_sequences_code)
        random_numbers_binary_huffman = n_huffman_compressor(random_numbers_binary[:repeat_number[n - 1]], n, all_n_p_sequences_code)
        lengths[index][n - 1] = int(len(random_numbers_binary_huffman) * 1000 / repeat_number[n - 1])
    index += 1
print(lengths)
fig = plt.figure(figsize=(15, 6))
fig.add_subplot(131)
plt.plot(block_lengths, lengths[0, :])
plt.plot(block_lengths, np.ones(len(block_lengths)) * 1000 * (-0.5 * np.log2(0.5) - 0.5 * np.log2(0.5)))
plt.ylim(950, 1050)
plt.xlabel('n')
plt.ylabel('length of compressed sequence')
plt.title('p = 0.5')
fig.add_subplot(132)
plt.plot(block_lengths, lengths[1, :])
plt.plot(block_lengths, np.ones(len(block_lengths)) * 1000 * (-0.75 * np.log2(0.75) - 0.25 * np.log2(0.25)))
plt.ylim(750, 1050)
plt.xlabel('n')
plt.ylabel('length of compressed sequence')
plt.title('p = 0.75')
fig.add_subplot(133)
plt.plot(block_lengths, lengths[2, :])
plt.plot(block_lengths, np.ones(len(block_lengths)) * 1000 * (-0.97* np.log2(0.97) - 0.03 * np.log2(0.03)))
plt.ylim(150, 1050)
plt.xlabel('n')
plt.ylabel('length of compressed sequence')
plt.title('p = 0.97')
plt.show()