from pymol import cmd

# modify this file to add your own filter

cmd.alter_state(1, 'all', 'x**=10; y**=10; z**=10')