pub mod math_eval;
pub mod menu;
pub mod newton;
pub mod placeholder;
pub mod theme;

use std::fmt::{self, Display, Formatter};

use strum::VariantArray;

pub const MODULES: &[Modules] = Modules::VARIANTS;

#[derive(Debug, Default, Clone, Copy, VariantArray)]
pub enum Modules {
    #[default]
    ExprEval,
    Newton,
    Placeholder,
}

impl Display for Modules {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Modules::ExprEval => write!(f, "Mathematical Evaluator"),
            Modules::Newton => write!(f, "Newton's Method"),
            Modules::Placeholder => write!(f, "Placeholder"),
        }
    }
}
