use ratatui::{
    style::{Color, Style},
    widgets::{Block, Borders},
};

pub struct Theme {}

impl Theme {
    pub const INSERT_FG: Style = Style::new().fg(Color::Green);
    pub const INSERT_BG: Style = Style::new().bg(Color::Green);
    pub const SELECTED_FG: Style = Style::new().fg(Color::Blue);
    pub const SELECTED_BG: Style = Style::new().bg(Color::Blue);
    pub const UNSELECTED_FG: Style = Style::new().fg(Color::Gray);
    pub const UNSELECTED_BG: Style = Style::new().bg(Color::DarkGray);

    pub fn block(style: Style) -> Block<'static> {
        Block::default().borders(Borders::ALL).border_style(style)
    }
}
