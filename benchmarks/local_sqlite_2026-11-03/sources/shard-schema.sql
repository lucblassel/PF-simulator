-- Schema setup and inserts share one transaction, including the version marker.
BEGIN IMMEDIATE;

PRAGMA user_version = 3;

CREATE TABLE pairs (
    id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, shard TEXT NOT NULL,
    tree_name TEXT NOT NULL, tree_offset INTEGER NOT NULL, tree_size INTEGER NOT NULL,
    sequence_name TEXT NOT NULL, sequence_offset INTEGER NOT NULL, sequence_size INTEGER NOT NULL
);

CREATE TABLE sample_metadata (
    name TEXT PRIMARY KEY REFERENCES pairs(name), candidate_id INTEGER NOT NULL,
    metadata_json TEXT NOT NULL
);

CREATE TABLE dataset_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
INSERT INTO dataset_metadata VALUES ('simulation_metadata_version', '1');
