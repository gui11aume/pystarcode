#!/usr/bin/en python
# -*- coding:utf-8 -*-


import sys
import gc
import Queue
import multiprocessing

from collections import defaultdict

class DupSeq(Exception):
   """Alert for sequence duplication."""
   pass

class MaxDistExceeded(Exception):
   """Allows to escape useless computation."""
   pass


def comp(seq1, seq2, maxdist, Npen):
   """'seq1' has to be longer than 'seq2'."""
   dist = 0
   for i in range(len(seq2)):
      if seq1[i] == 'N' or seq2[i] == 'N':
         dist += Npen
      elif seq1[i] == seq2[i]:
         continue
      else:
         dist += 1
      if dist > maxdist:
         raise MaxDistExceeded
   return dist


def ppref(string1, string2):
   """Return the common prefix between two strings, and the two
   corresponding suffixes. 'string1' has to be the longest."""

   for num,a in enumerate(string2):
      if string1[num] != a: break
   return (string1[:num], string1[num:], string2[num:])



class seqTrieNode(object):
   """Nodes of a seqTrie have 3 attributes:
       'subseq'  : the characters on a seq path
       'data'    : replacement characters ('None' if node is not tail)
       'children': list of 'seqTrieNode' children objects."""

   def __init__(self, parent, subseq='', data=None):
      """seqTrie nodes are instantiated with no child."""
      self.parent = parent
      self.children = []
      self.subseq = subseq
      self.data = data

   @property
   def ancestry(self):
      """Orgin first ancestry iterator."""
      if self.parent is not None:
         for ancestor in self.parent.ancestry:
            yield ancestor
      yield self

   @property
   def seq(self):
      if self.data is None:
         return ''
      return ''.join([node.subseq for node in self.ancestry])

   def fork(self, seq, data):
      """A method for readability. Use 'fork' to add a node that
      shares a prefix with current node to the trie."""

      (pref, suff_1, suff_2) = ppref(self.subseq, seq)
      if pref == '':
         leaf = seqTrieNode(parent=self.parent, subseq=seq, data=data)
         self.parent.children.append(leaf)
      else:
         internal = seqTrieNode(parent=self.parent, subseq=pref)
         leaf = seqTrieNode(parent=internal, subseq=suff_2, data=data)
         internal.children = [self, leaf]
         self.subseq = suff_1
         self.parent.children.append(internal)
         self.parent.children.remove(self)
         self.parent = internal
      return leaf


class seqTrie(object):
   """A radix trie is a collection of nodes linked by parent-child
   relationships. This class contais the methods to build and search
   the trie."""

   def __init__(self):
      """Radix tries are instantiated with an empty 'root' node."""
      self.root = seqTrieNode(parent=None)

   def walk(self, node=None, depth=-1):
      if node is None: node = self.root
      yield depth,node
      for child in node.children:
         for newdepth,newnode in self.walk(child,depth+1):
            yield newdepth,newnode

   def search(self, seq, node=None, maxdist=0, Npen=.251):
      """Depth-first iterator that returns all the matches to a
      query sequence in the seqTrie with a distance less than or
      equal to the specified max."""

      if node is None: node = self.root
      for child in node.children:
         try:
            dist = comp(seq, child.subseq, maxdist, Npen)
         except MaxDistExceeded:
            # Max distance exceeded: no need to go deeper.
            continue
         if child.data is not None:
            yield (child.seq, child.data)
         subseq = seq[len(child.subseq):]
         for hit in self.search(subseq, child, maxdist-dist):
            yield hit
      return

   def insert_sorted(self, item, node):
      """Insert the item to the trie (works only for sorted items)."""
      seq, data = item
      for ancestor in node.ancestry:
         if seq.startswith(ancestor.subseq):
            seq = seq[len(ancestor.subseq):]
            continue
         else:
            break
      return ancestor.fork(seq, data)

   @classmethod
   def build_from_dict(self, pairdict):
      """Conveniently build a radix trie from key/data items."""
      # The list 'append' method has a bug in Python.
      # See http://bugs.python.org/issue4074
      # Disabling garbage collection while doing multiple 'append'
      # of complex objects can increase performance.
      gc_was_initially_enabled = gc.isenabled()
      try:
         gc.disable()
         trie = seqTrie()
         # Check that no key is duplicated.
         if len(pairdict.keys()) != len(set(pairdict.keys())):
            raise KeyDuplicated()
         node = trie.root
         sorted_items = sorted(pairdict.items())
         seq,data = sorted_items.pop(0)
         node = seqTrieNode(trie.root, seq, data)
         trie.root.children.append(node)
         for item in sorted_items:
            node = trie.insert_sorted(item, node)
      finally:
         if gc_was_initially_enabled: gc.enable()
      return trie

class searchThread(multiprocessing.Process):
   def __init__(self, trie, queuelock, nodelock, queue, nodes, **args):
      self.trie = trie
      self.queuelock = queuelock
      self.nodelock = nodelock
      self.queue = queue
      self.nodes = nodes
      self.args = args
      multiprocessing.Process.__init__(self)

   def run(self):
      """Pure side effect fuction that updates the 'nodes' dictionary
      in place."""

      while True:
         self.queuelock.acquire()
         try:
            item = queue.get_nowait()
         except Queue.Empty:
            return
         finally:
            self.queuelock.release()

         left = '%s:%d' % item
         self.nodelock.acquire()
         self.nodes[left] = []
         self.nodelock.release()
         for hit in self.trie.search(item[0], **self.args):
            right = '%s:%d' % hit
            # Prevent cycles in the graph.
            if left < right:
               self.nodelock.acquire()
               self.nodes[left] += [right]
               self.nodelock.release()
         self.queue.task_done()


if __name__ == '__main__':
   sequences = defaultdict(int)
   with open(sys.argv[1]) as f:
      sys.stderr.write('reading file\n')
      for line in f:
         sequences[line.rstrip()] += 1

   sys.stderr.write('building trie\n')
   trie = seqTrie.build_from_dict(sequences)

   sys.stderr.write('matching pairs\n')

   queuelock = multiprocessing.Lock()
   nodelock = multiprocessing.Lock()
   queue = multiprocessing.JoinableQueue(-1)
   manager = multiprocessing.Manager()
   nodes_ = manager.dict()
   
   for item in sequences.items():
      queue.put_nowait(item)

   processes = [
      searchThread(
         trie = trie,
         queuelock = queuelock,
         nodelock = nodelock,
         queue = queue,
         nodes = nodes_,
         maxdist = 2,
         Npen = .251
      ) for i in range(6)
   ]

   for process in processes:
      process.start()
      sys.stderr.write('starting process %d\n' % process.pid)
   queue.join()

   # The class 'manager.dict' shares only part of the 'dict'
   # interface, which causes the following operations to crash.
   nodes = {}
   for node in nodes_.keys():
      nodes[node] = nodes_[node]
   del(nodes_)

   sys.stderr.write('writing pair file\n')
   with open('pairs.txt', 'w') as f:
      for source,targets in nodes.items():
         for target in targets:
            f.write('%s\t%s\n' % (source,target))

   sys.stderr.write('writing seq file\n')
   with open('seq.txt', 'w') as f:
      for node in nodes.keys():
         f.write(node + '\n')

   sys.stderr.write('finding connected components\n')
   def getCC(node, graph, S):
      S.update([node])
      for newnode in graph.pop(node, []):
         getCC(newnode, graph, S)

   CC = []
   while nodes:
      S = set()
      node = nodes.iterkeys().next()
      getCC(node, nodes, S)
      CC.append(S)

   for cc in CC:
      sys.stdout.write(', '.join(cc) + '\n')

