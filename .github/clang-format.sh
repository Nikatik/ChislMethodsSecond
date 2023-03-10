#!/usr/bin/env bash
clang-format --style="{BasedOnStyle: llvm, IndentWidth: 4, AllowShortFunctionsOnASingleLine: Empty, AlignConsecutiveAssignments: Consecutive, AlignArrayOfStructures: Right, AlignTrailingComments: true, AllowShortBlocksOnASingleLine: Empty, AllowShortIfStatementsOnASingleLine: AllIfsAndElse, IndentCaseLabels: true, IndentCaseBlocks: true, IndentAccessModifiers: true, BinPackArguments: false, BinPackParameters: false, BitFieldColonSpacing: After, BreakBeforeBraces: Allman, EmptyLineBeforeAccessModifier: Always, PointerAlignment: Left, ReferenceAlignment: Left, SeparateDefinitionBlocks: Always, SpaceBeforeParens: NonEmptyParentheses, SpacesBeforeTrailingComments: 8, Standard: Latest}" "$@"
