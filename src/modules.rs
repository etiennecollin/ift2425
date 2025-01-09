pub mod math_eval;
pub mod menu;
pub mod placeholder;
pub mod theme;

use std::fmt::{self, Display, Formatter};

use strum::VariantArray;

pub const MODULES: &[Modules] = Modules::VARIANTS;

#[derive(Debug, Default, Clone, Copy, VariantArray)]
pub enum Modules {
    #[default]
    ExprEval,
    Placeholder,
}

impl Display for Modules {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Modules::ExprEval => write!(f, "Expression Evaluator"),
            Modules::Placeholder => write!(f, "Placeholder"),
        }
    }
}
