use ratatui::{
    buffer::Buffer,
    layout::{Constraint, Layout, Rect},
    style::Stylize,
    widgets::{Paragraph, Widget},
};

use super::theme::Theme;

#[derive(Debug, Default)]
pub struct PlaceholderModule {}

impl PlaceholderModule {}
impl Widget for &PlaceholderModule {
    fn render(self, area: Rect, buf: &mut Buffer) {
        // Create the title area
        let vertical = Layout::vertical([
            Constraint::Fill(1),
            Constraint::Fill(1),
            Constraint::Length(3),
            Constraint::Length(3),
        ]);

        let [expr_area, x_area, result_area, error_area] = vertical.areas(area);

        let result_block = Theme::block(Theme::UNSELECTED_FG).title("Placeholder Title");
        let result = Paragraph::new("Placeholder text")
            .block(result_block)
            .centered()
            .bold();

        let status_block = Theme::block(Theme::UNSELECTED_FG).title("PLaceholder Title");
        let status = Paragraph::new("Placeholder text")
            .block(status_block)
            .centered()
            .bold();

        let expr_block = Theme::block(Theme::UNSELECTED_FG).title("Placeholder Title");
        let expr = Paragraph::new("Placeholder text")
            .block(expr_block)
            .centered()
            .bold();

        let x_block = Theme::block(Theme::UNSELECTED_FG).title("Placeholder Title");
        let x = Paragraph::new("Placeholder text")
            .block(x_block)
            .centered()
            .bold();

        // Render everything
        expr.render(expr_area, buf);
        x.render(x_area, buf);
        result.render(result_area, buf);
        status.render(error_area, buf);
    }
}
