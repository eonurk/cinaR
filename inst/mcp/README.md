# cinaR MCP starter spec

This folder contains starter MCP tool metadata for `cinaR`.

- `cinaR-mcp-tools.json`: tool definitions with input/output JSON Schemas.

Notes:

- This file is a schema/spec template, not a running MCP server by itself.
- Implementations should map file paths/arguments in the schemas to calls of:
  - `cinaR()`
  - `prep_scATAC_cinaR()`
  - `prep_scATAC_seurat()`
