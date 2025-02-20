use std::fmt::{self, Display, Formatter};

use crossterm::event::KeyEvent;
use ratatui::{
    buffer::Buffer,
    layout::{Constraint, Layout, Rect},
    style::Stylize,
    widgets::{Paragraph, Widget},
};
use tui_textarea::TextArea;
use utils::newton::newton;

use super::theme::Theme;

const EXPR_AREA_TITLE: &str = "X";
const CONTEXT_AREA_TITLE: &str = "Array";

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
pub struct NewtonModule {
    expr_area: TextArea<'static>,
    context_area: TextArea<'static>,
    result: f64,
    compute_state: ComputeState,
    current_field: CurrentField,
}

impl Default for NewtonModule {
    fn default() -> Self {
        let mut expr_area = TextArea::default();
        expr_area.set_placeholder_text("Enter a valid array (e.g. 1, 2, 3, 10, -2.3)");
        expr_area.set_cursor_style(Theme::UNSELECTED_BG);
        expr_area.set_block(Theme::block(Theme::SELECTED_FG).title(EXPR_AREA_TITLE));

        let mut x_area = TextArea::default();
        x_area.set_placeholder_text("Enter a valid x (e.g. 3.26)");
        x_area.set_cursor_style(Theme::UNSELECTED_BG);
        x_area.set_block(Theme::block(Theme::UNSELECTED_FG).title(CONTEXT_AREA_TITLE));

        Self {
            expr_area,
            context_area: x_area,
            result: 0.0,
            compute_state: ComputeState::Idle,
            current_field: CurrentField::Expression,
        }
    }
}

impl NewtonModule {
    pub fn run(&mut self) {
        self.compute_state = ComputeState::Computing;
        match self.parse_expr() {
            Ok(()) => (),
            Err(()) => return,
        };
        self.compute_state = ComputeState::Done;
    }

    fn parse_expr(&mut self) -> Result<(), ()> {
        // Get the contents of the expression textarea
        let x_str = match self.expr_area.lines().first() {
            Some(line) => line,
            None => {
                self.compute_state = ComputeState::Error("Empty x".to_string());
                return Err(());
            }
        };

        // Parse the list of numbers into a Vec<f64>
        let x: f64 = match x_str.trim().parse::<f64>() {
            Ok(x) => x,
            Err(e) => {
                self.compute_state = ComputeState::Error(format!("x error: {}", e));
                return Err(());
            }
        };

        // Get the contents of the context textarea
        let array_str = match self.context_area.lines().first() {
            Some(line) => line,
            None => {
                self.compute_state = ComputeState::Error("Empty array".to_string());
                return Err(());
            }
        };

        // Parse the list of numbers into a Vec<f64>
        let array: Vec<f64> = match array_str
            .split(',')
            .map(|x| x.trim().parse::<f64>())
            .collect::<Result<Vec<f64>, _>>()
        {
            Ok(x) => x,
            Err(e) => {
                self.compute_state = ComputeState::Error(format!("Array error: {}", e));
                return Err(());
            }
        };

        self.result = newton(x, &array);
        Ok(())
    }

    pub fn input(&mut self, c: KeyEvent) {
        match self.current_field {
            CurrentField::Expression => self.expr_area.input(c),
            CurrentField::Value => self.context_area.input(c),
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
                    .set_block(Theme::block(Theme::UNSELECTED_FG).title(EXPR_AREA_TITLE));
            }
            CurrentField::Value => {
                self.context_area.set_cursor_style(Theme::UNSELECTED_BG);
                self.context_area
                    .set_block(Theme::block(Theme::UNSELECTED_FG).title(CONTEXT_AREA_TITLE));
            }
        }
    }

    fn field_select(&mut self) {
        match self.current_field {
            CurrentField::Expression => {
                self.expr_area.set_cursor_style(Theme::SELECTED_BG);
                self.expr_area
                    .set_block(Theme::block(Theme::SELECTED_FG).title(EXPR_AREA_TITLE));
            }
            CurrentField::Value => {
                self.context_area.set_cursor_style(Theme::SELECTED_BG);
                self.context_area
                    .set_block(Theme::block(Theme::SELECTED_FG).title(CONTEXT_AREA_TITLE))
            }
        }
    }

    pub fn insert_mode_enter(&mut self) {
        match self.current_field {
            CurrentField::Expression => {
                self.expr_area.set_cursor_style(Theme::INSERT_BG);
                self.expr_area
                    .set_block(Theme::block(Theme::INSERT_FG).title(EXPR_AREA_TITLE))
            }
            CurrentField::Value => {
                self.context_area.set_cursor_style(Theme::INSERT_BG);
                self.context_area
                    .set_block(Theme::block(Theme::INSERT_FG).title(CONTEXT_AREA_TITLE))
            }
        }
    }

    pub fn insert_mode_exit(&mut self) {
        self.field_select();
    }
}

impl Widget for &NewtonModule {
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
        self.context_area.render(x_area, buf);
        result.render(result_area, buf);
        status.render(error_area, buf);
    }
}
