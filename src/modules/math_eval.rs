use std::fmt::{self, Display, Formatter};

use crossterm::event::KeyEvent;
use ratatui::{
    buffer::Buffer,
    layout::{Constraint, Layout, Rect},
    style::Stylize,
    widgets::{Paragraph, Widget},
};
use tui_textarea::TextArea;

use super::theme::Theme;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
enum ComputeState {
    #[default]
    Idle,
    Computing,
    Done,
    Error(String),
}

impl Display for ComputeState {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ComputeState::Idle => write!(f, "Idle"),
            ComputeState::Computing => write!(f, "Computing..."),
            ComputeState::Done => write!(f, "Done"),
            ComputeState::Error(e) => write!(f, "Error: {}", e),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum CurrentField {
    #[default]
    Expression,
    Value,
}

#[derive(Debug)]
pub struct MathEvalModule {
    expr_area: TextArea<'static>,
    x_area: TextArea<'static>,
    result: f64,
    compute_state: ComputeState,
    current_field: CurrentField,
}

impl Default for MathEvalModule {
    fn default() -> Self {
        let mut expr_area = TextArea::default();
        expr_area.set_placeholder_text("Enter a valid expression (e.g. pi * x^e + 1)");
        expr_area.set_cursor_style(Theme::UNSELECTED_BG);
        expr_area.set_block(Theme::block(Theme::SELECTED_FG).title("Mathematical Expression"));

        let mut x_area = TextArea::default();
        x_area.set_placeholder_text(
            "Enter a valid expression that evaluates to a number (e.g. 1.56 * ln(pi)^asin(e/3))",
        );
        x_area.set_cursor_style(Theme::UNSELECTED_BG);
        x_area.set_block(Theme::block(Theme::UNSELECTED_FG).title("Value of X"));

        Self {
            expr_area,
            x_area,
            result: 0.0,
            compute_state: ComputeState::Idle,
            current_field: CurrentField::Expression,
        }
    }
}

impl MathEvalModule {
    pub fn run(&mut self) {
        let (func, x) = match self.parse_expr() {
            Ok((func, x)) => (func, x),
            Err(()) => return,
        };

        self.compute_state = ComputeState::Computing;
        self.result = func(x);
        self.compute_state = ComputeState::Done;
    }

    fn parse_expr(&mut self) -> Result<(impl Fn(f64) -> f64, f64), ()> {
        // Get the first line of the textarea
        let line = match self.expr_area.lines().first() {
            Some(line) => line,
            None => {
                self.compute_state = ComputeState::Error("Empty expression".to_string());
                return Err(());
            }
        };

        // Parse the expression into a mathematical expression
        let expr: meval::Expr = match line.parse() {
            Ok(expr) => expr,
            Err(e) => {
                self.compute_state = ComputeState::Error(format!("Expression error: {}", e));
                return Err(());
            }
        };

        // Use the character `x` as the variable for the function expression
        let func = match expr.bind("x") {
            Ok(func) => func,
            Err(e) => {
                self.compute_state = ComputeState::Error(format!("Expression error: {}", e));
                return Err(());
            }
        };

        // Get the value of `x` from the module
        let line = match self.x_area.lines().first() {
            Some(line) => line,
            None => {
                self.compute_state = ComputeState::Error("Empty x value".to_string());
                return Err(());
            }
        };

        let x: f64 = match meval::eval_str(line) {
            Ok(x) => x,
            Err(e) => {
                self.compute_state = ComputeState::Error(format!("X value error: {}", e));
                return Err(());
            }
        };

        Ok((func, x))
    }

    pub fn input(&mut self, c: KeyEvent) {
        match self.current_field {
            CurrentField::Expression => self.expr_area.input(c),
            CurrentField::Value => self.x_area.input(c),
        };
    }

    pub fn next_field(&mut self) {
        self.field_unselect();
        match self.current_field {
            CurrentField::Expression => self.current_field = CurrentField::Value,
            CurrentField::Value => self.current_field = CurrentField::Expression,
        }
        self.field_select();
    }

    pub fn prev_field(&mut self) {
        self.field_unselect();
        match self.current_field {
            CurrentField::Expression => self.current_field = CurrentField::Value,
            CurrentField::Value => self.current_field = CurrentField::Expression,
        }
        self.field_select();
    }

    fn field_unselect(&mut self) {
        match self.current_field {
            CurrentField::Expression => {
                self.expr_area.set_cursor_style(Theme::UNSELECTED_BG);
                self.expr_area
                    .set_block(Theme::block(Theme::UNSELECTED_FG).title("Mathematical Expression"));
            }
            CurrentField::Value => {
                self.x_area.set_cursor_style(Theme::UNSELECTED_BG);
                self.x_area
                    .set_block(Theme::block(Theme::UNSELECTED_FG).title("Value of X"))
            }
        }
    }

    fn field_select(&mut self) {
        match self.current_field {
            CurrentField::Expression => {
                self.expr_area.set_cursor_style(Theme::SELECTED_BG);
                self.expr_area
                    .set_block(Theme::block(Theme::SELECTED_FG).title("Mathematical Expression"));
            }
            CurrentField::Value => {
                self.x_area.set_cursor_style(Theme::SELECTED_BG);
                self.x_area
                    .set_block(Theme::block(Theme::SELECTED_FG).title("Value of X"))
            }
        }
    }

    pub fn insert_mode_enter(&mut self) {
        match self.current_field {
            CurrentField::Expression => {
                self.expr_area.set_cursor_style(Theme::INSERT_BG);
                self.expr_area
                    .set_block(Theme::block(Theme::INSERT_FG).title("Mathematical Expression"));
            }
            CurrentField::Value => {
                self.x_area.set_cursor_style(Theme::INSERT_BG);
                self.x_area
                    .set_block(Theme::block(Theme::INSERT_FG).title("Value of X"))
            }
        }
    }

    pub fn insert_mode_exit(&mut self) {
        self.field_select();
    }
}

impl Widget for &MathEvalModule {
    fn render(self, area: Rect, buf: &mut Buffer) {
        // Create the title area
        let vertical = Layout::vertical([
            Constraint::Fill(1),
            Constraint::Fill(1),
            Constraint::Length(3),
            Constraint::Length(3),
        ]);

        let [expr_area, x_area, result_area, error_area] = vertical.areas(area);

        let result_block = Theme::block(Theme::UNSELECTED_FG).title("Result");
        let result = Paragraph::new(format!("{}", self.result))
            .block(result_block)
            .centered()
            .bold();

        let status_block = Theme::block(Theme::UNSELECTED_FG).title("Status");
        let status = Paragraph::new(format!("{}", self.compute_state))
            .block(status_block)
            .centered()
            .bold();

        // Render everything
        self.expr_area.render(expr_area, buf);
        self.x_area.render(x_area, buf);
        result.render(result_area, buf);
        status.render(error_area, buf);
    }
}
