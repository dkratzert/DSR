#/usr/bin/env python
#-*- encoding: utf-8 -*-
#möp
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import os
from collections import OrderedDict

import misc

feedpath = os.path.abspath("setup/Output/changelog.txt")
#outpath = 'W:/htdocs/data/feed.rss'
outpath = './feed.rss'

head = r"""
<?xml version="1.0" encoding="utf-8"?>

<rss version="2.0">

  <channel>
    <title>DSR Changelog</title>
    <link>https://www.xs3.uni-freiburg.de/research/dsr</link>
    <description>Recent changes made to DSR development.</description>
    <language>en-en</language>
    <copyright>Daniel Kratzert</copyright>
    <pubDate>05.01.2017</pubDate>
"""
item = r"""
    <item>
      <title> {} </title>
      <description>{}</description>
      <link>https://www.xs3.uni-freiburg.de/research/dsr</link>
    </item>
"""

foot = r"""
  </channel>

</rss>
"""

def readfile(feedpath):
    with open(feedpath, 'r') as f:
        #changelst = f.readlines()
        changelst = f.read().splitlines()
    return changelst

def versionsearcher(lines):
    """
    Searches for versions in changelog
    """
    occourences = misc.find_multi_lines(lines, '^\-[0-9]{0,5}')
    return occourences


def get_textblocks(lines, occ):
    """
    parses textblocks from changelog
    """
    itemlist = OrderedDict()
    for num, i in enumerate(occ, 1):
        try:
            rawitem = ' '.join([l.strip() for l in lines[i:occ[num]]])
            #print(rawitem)
        except IndexError:
            #print('index error')
            rawitem = ' '.join([l.strip() for l in lines[i:]])
        if rawitem.startswith('-'):
            rawitem = rawitem.replace('\t', ' ')
            line = rawitem.split(' ')
            version = line[0].strip('-')
            if num == 1:
                print('Newest version:', version)
            entry = ' '.join(line[1:]) # text
            #print(version, entry)
            itemlist[version] = entry
    return itemlist

def xml_creator(items):
    """
    put all together for rss output
    """
    rss = " "
    rss = rss+head
    for it in items:
        rss = rss+item.format(it, items[it])
    rss = rss+foot
    return rss

def filewriter(rss):
    with open(os.path.abspath(outpath), 'w') as f:
        f.write(rss)


if __name__ == "__main__":
    filedata = readfile(feedpath)
    occ = versionsearcher(filedata)
    items = get_textblocks(filedata, occ)
    rss = xml_creator(items)
    #print(rss)
    filewriter(rss)
    print("new feed written")