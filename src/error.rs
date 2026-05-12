use std::path::PathBuf;

/// Errors that can occur in rustar-aligner.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid parameter: {0}")]
    Parameter(String),

    #[error("I/O error: {source} ({path})")]
    Io {
        source: std::io::Error,
        path: PathBuf,
    },

    #[error("FASTA parsing error: {0}")]
    Fasta(String),

    #[error("genome index error: {0}")]
    Index(String),

    #[error("alignment error: {0}")]
    Alignment(String),

    #[error("GTF parsing error: {0}")]
    Gtf(String),

    #[error("chimeric detection error: {0}")]
    Chimeric(String),
}

impl Error {
    /// Convenience for wrapping an `io::Error` with a path context.
    pub fn io(source: std::io::Error, path: impl Into<PathBuf>) -> Self {
        Self::Io {
            source,
            path: path.into(),
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(err: std::io::Error) -> Self {
        Self::Io {
            source: err,
            path: PathBuf::from("<unknown>"),
        }
    }
}
