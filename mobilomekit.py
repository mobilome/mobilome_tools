#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   mobilomekit.py
@Time    :   2025/01/05 21:29:55
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com
'''

# here put the import lib
import argparse
from mobilomekit import script1, script2, script3

def main():
    parser = argparse.ArgumentParser(description="mobilomekit - A toolkit for the mobilome laboratory.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subcommand for script1
    parser_script1 = subparsers.add_parser('script1', help='Run script1', description='This command runs script1, which performs a specific analysis on the input data.')
    parser_script1.add_argument('--input', required=True, help='Input file for script1')
    parser_script1.add_argument('--output', help='Output file for script1')

    # Subcommand for script2
    parser_script2 = subparsers.add_parser('script2', help='Run script2')
    parser_script2.add_argument('--input', required=True, help='Input file for script2')
    parser_script2.add_argument('--output', help='Output file for script2')

    # Subcommand for script3
    parser_script3 = subparsers.add_parser('script3', help='Run script3')
    parser_script3.add_argument('--input', required=True, help='Input file for script3')
    parser_script3.add_argument('--output', help='Output file for script3')

    args = parser.parse_args()

    if args.command == 'script1':
        script1.run(args.input, args.output)
    elif args.command == 'script2':
        script2.run(args.input, args.output)
    elif args.command == 'script3':
        script3.run(args.input, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()