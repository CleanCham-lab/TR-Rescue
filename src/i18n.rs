// TR-Rescue — Copyright (c) 2026 Seung Won Chae. SPDX-License-Identifier: MIT
//! Korean / English text pairs. The report carries both; TSV files and the CLI use `--lang`.

use serde::Serialize;

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Tr {
    pub ko: String,
    pub en: String,
}

impl Tr {
    pub fn new(ko: impl Into<String>, en: impl Into<String>) -> Tr {
        Tr { ko: ko.into(), en: en.into() }
    }
    pub fn get(&self, lang: &str) -> &str {
        if lang == "en" { &self.en } else { &self.ko }
    }
    pub fn is_empty(&self) -> bool {
        self.ko.is_empty() && self.en.is_empty()
    }
}

/// Pick the Korean or English string.
pub fn t<'a>(lang: &str, ko: &'a str, en: &'a str) -> &'a str {
    if lang == "en" { en } else { ko }
}

/// Join a list of pairs with separators per language.
pub fn join(v: &[Tr], sep: &str) -> Tr {
    Tr::new(v.iter().map(|x| x.ko.as_str()).collect::<Vec<_>>().join(sep), v.iter().map(|x| x.en.as_str()).collect::<Vec<_>>().join(sep))
}
