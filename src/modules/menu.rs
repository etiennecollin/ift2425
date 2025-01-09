use std::cell::RefCell;

use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Style, Stylize},
    widgets::{HighlightSpacing, List, ListDirection, ListState, StatefulWidget, Widget},
};
use strum::VariantArray;

use super::{theme::Theme, Modules};

#[derive(Debug)]
pub struct MenuModule {
    items: Vec<String>,
    state: RefCell<ListState>,
    current_index: usize,
}

impl Default for MenuModule {
    fn default() -> Self {
        // Create a menu item for each module
        let items: Vec<String> = Modules::VARIANTS
            .iter()
            .map(|module| module.to_string())
            .collect();

        // Set the first module as selected
        let mut state = ListState::default();
        state.select(Some(0));

        Self {
            items,
            state: state.into(),
            current_index: 0,
        }
    }
}

impl MenuModule {
    pub fn next_module(&mut self) {
        // Wrap around to the first module if we reach the end
        if self.current_index == self.items.len() - 1 {
            self.state.borrow_mut().select(Some(0));
            self.current_index = 0;
        } else {
            self.state.borrow_mut().scroll_down_by(1);
            self.current_index += 1;
        }
    }

    pub fn prev_module(&mut self) {
        if self.current_index == 0 {
            self.state.borrow_mut().select(Some(self.items.len() - 1));
            self.current_index = self.items.len() - 1;
        } else {
            self.state.borrow_mut().scroll_up_by(1);
            self.current_index -= 1;
        }
    }

    pub fn get_selected_module(&self) -> usize {
        self.state
            .borrow()
            .selected()
            .expect("A module should always be selected. This should never happen.")
    }
}

impl Widget for &MenuModule {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let block = Theme::block(Theme::UNSELECTED_FG).title("Menu");
        let list = List::new(self.items.clone())
            .block(block)
            .highlight_spacing(HighlightSpacing::Always)
            .highlight_style(Style::new().reversed())
            .highlight_symbol(">>")
            .direction(ListDirection::TopToBottom);

        StatefulWidget::render(list, area, buf, &mut self.state.borrow_mut());
    }
}
