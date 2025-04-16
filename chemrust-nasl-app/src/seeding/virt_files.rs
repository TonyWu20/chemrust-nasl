use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct FileModel {
    path: PathBuf,
    data: Vec<u8>,
}

impl FileModel {
    pub fn new(path: PathBuf, data: Vec<u8>) -> Self {
        Self { path, data }
    }

    pub fn path(&self) -> &PathBuf {
        &self.path
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }
}

#[derive(Debug, Clone)]
pub struct ContentStorage {
    files: Vec<FileModel>,
}

impl ContentStorage {
    pub fn new(files: Vec<FileModel>) -> Self {
        Self { files }
    }
    pub fn add_file(&mut self, file: FileModel) {
        self.files.push(file);
    }

    pub fn files(&self) -> &[FileModel] {
        &self.files
    }
}
